from collections import namedtuple
import logging
from pathlib import Path
import hashlib
import yaml
from tempfile import TemporaryDirectory
import tarfile

from django.utils import timezone

from django.db import IntegrityError
from django.db import transaction


from django.contrib.auth import get_user_model


from viewer.models import (
    Compound,
    Project,
    Xtalform,
    XtalformSite,
    QuatAssembly,
    CanonSite,
    CanonSiteConf,
    SiteObservation,
    Experiment,
    ExperimentUpload,
)

from viewer.target_set_upload import get_create_target

from django.conf import settings

logger = logging.getLogger(__name__)

METADATA_FILE = "meta_aligner.yaml"
TARGET_LOADER_DATA = "target_loader_data"

# process all metadata yaml file blocks into following structure:
# { some_id: MetadataObjects, ...}
# reason being, quite often I need to refer to these by
# some alternative id. with data struct in namedtuple, I'm able
# to create temporary dicts with key that's needed atm
MetadataObjects = namedtuple("MetadataObjects", "instance data")


class TargetLoader:
    def __init__(self, data_bundle: str, tempdir: str):
        self.data_bundle = data_bundle
        self.bundle_name = Path(data_bundle).stem
        self.target_name = self.bundle_name.split("-")[0]
        self.bundle_path = Path(settings.MEDIA_ROOT, data_bundle)
        self.tempdir = tempdir
        self.raw_data = Path(self.tempdir).joinpath(self.bundle_name)

        self._target_root = self.raw_data.joinpath(self.target_name).joinpath(
            "upload_1"
        )

    @property
    def target_root(self):
        return self._target_root

    def _load_yaml_blocks(self, yaml_file, blocks) -> list:
        try:
            with open(yaml_file, "r", encoding="utf-8") as file:
                meta_aligner = yaml.safe_load(file)
        except FileNotFoundError as exc:
            msg = f"{METADATA_FILE} file not found in data archive!"
            logger.error(msg)
            raise FileNotFoundError(msg) from exc

        validation_errors = []
        error_text = "'{}' section missing in input file"
        result = []
        for block in blocks:
            try:
                result.append(meta_aligner[block])
            except KeyError:
                msg = error_text.format(block)
                validation_errors.append(msg)
                logger.error(msg)

        # wonder if it's worth throwing a custom exception here.. easier
        # to scan the logs. then again, if errors are passed to user,
        # inspecting logs isn't really necessary?
        if validation_errors:
            raise KeyError(validation_errors)

        return result

    def _process_experiment(self, experiment_upload, protein_name, data):
        """Create Experiment model instance from data.

        Incoming data format:
        type: <type>
        last_updated: <datetime>
        crystallographic_files:
          xtal_pdb: {file: <file path>, sha256: <hash>}
          xtal_mtz: {file: <file path>, sha256: <hash>}
          ligand_cif: {file: <file path>, sha256: <hash>, smiles: <smiles>}
          panddas_event_files:
          - {file: <file path>, sha256: <hash>,
            model: <int>, chain: <char[1]>, res: <int>, index: <int>, bdc: <float>}
          - {file: <file path>, sha256: <hash>,
            model: <int>, chain: <char[1]>, res: <int>, index: <int>, bdc: <float>}
        status: <status>
        assigned_xtalform: <int>
        aligned_files: [...]

        Manages to save all references to other tables
        """
        logger.debug("Creating experiment object: %s", protein_name)

        files = self._check_file_struct(
            self.target_root, data["crystallographic_files"]
        )

        experiment = Experiment(
            experiment_upload=experiment_upload,
            code=protein_name,
            status=1,
            version=1,
            type=1 if data["type"] is "manual" else 0,  # FIXME
            pdb_info=files["xtal_pdb"],
            mtz_info=files["xtal_mtz"],
            cif_info=files.get("cif_info", None),  # this may be missing from the struct
            # this doesn't seem to be present
            # pdb_sha256=
        )
        try:
            experiment.save()
        except IntegrityError as exc:
            msg = f"Failed to save Experiment: {protein_name}"
            logger.error(msg)
            raise IntegrityError(msg) from exc

        return MetadataObjects(
            instance=experiment, data=data.get("aligned_files", None)
        )

    def _process_compound(self, protein_name, data):
        """Create Compound model instance from data.

        Incoming data format:
        xtal_pdb: {file: <file path>, sha256: <hash>}
        xtal_mtz: {file: <file path>, sha256: <hash>}
        ligand_cif: {file: <file path>, sha256: <hash>, smiles: <smiles>}
        panddas_event_files:
        - {file: <file path>, sha256: <hash>,
          model: <int>, chain: <char[1]>, res: <int>, index: <int>, bdc: <float>}
        - {file: <file path>, sha256: <hash>,
          model: <int>, chain: <char[1]>, res: <int>, index: <int>, bdc: <float>}
        """

        # TODO: there's a method, calc_cpd, could I use that? as I
        # understand, they don't need to be unique anymore, so I can
        # just remove the uniqueness validation?

        logger.debug("Creating Compound object: %s", protein_name)
        compound = Compound(
            smiles=data["ligand_cif"]["smiles"],
            # mol_log_p = Chem.Crippen.MolLogP(m),
            # mol_wt = float(Chem.rdMolDescriptors.CalcExactMolWt(m)),
            # heavy_atom_count = Chem.Lipinski.HeavyAtomCount(m),
            # heavy_atom_mol_wt = float(Descriptors.HeavyAtomMolWt(m)),
            # nhoh_count = Chem.Lipinski.NHOHCount(m),
            # no_count = Chem.Lipinski.NOCount(m),
            # num_h_acceptors = Chem.Lipinski.NumHAcceptors(m),
            # num_h_donors = Chem.Lipinski.NumHDonors(m),
            # num_het_atoms = Chem.Lipinski.NumHeteroatoms(m),
            # num_rot_bonds = Chem.Lipinski.NumRotatableBonds(m),
            # num_val_electrons = Descriptors.NumValenceElectrons(m),
            # ring_count = Chem.Lipinski.RingCount(m),
            # tpsa = Chem.rdMolDescriptors.CalcTPSA(m),
            mol_log_p=0.0,
            mol_wt=0.0,
            heavy_atom_count=0,
            heavy_atom_mol_wt=0,
            nhoh_count=0,
            no_count=0,
            num_h_acceptors=0,
            num_h_donors=0,
            num_het_atoms=0,
            num_rot_bonds=0,
            num_val_electrons=0,
            ring_count=0,
            tpsa=0.0,
        )
        try:
            compound.save()
        except IntegrityError as exc:
            msg = f"Failed to save Compound: {protein_name}"
            logger.error(msg)
            raise IntegrityError(msg) from exc
        except KeyError as exc:
            # this means ligand info missing
            raise KeyError from exc

        # data basically just contains file locations, need to copy them
        # later
        return MetadataObjects(instance=compound, data=data)

    def _process_xtalform(self, experiment, quat_assembly, idx, data):
        """Create Xtalform model instance from data.

        Incoming data format:
        <idx>:
          xtalform_ref: <ref>
          xtalform_space_group: <space group>
          xtalform_cell: <cell info>

        Saves all references to other tables (QuatAssembly and Experiment).
        """
        logger.debug("Creating Xtalform object: %s", data["xtalform_ref"])
        xtalform = Xtalform(
            experiment=experiment,
            quat_assembly=quat_assembly,
            name=data["xtalform_ref"],
            space_group=data["xtalform_space_group"],
            unit_cell_info=data["xtalform_cell"],
            xtalform_id=idx,
        )
        try:
            xtalform.save()
        except IntegrityError as exc:
            msg = f"Failed to save Xtalform: {data['xtalform_ref']}"
            logger.error(msg)
            raise IntegrityError(msg) from exc

        return MetadataObjects(instance=xtalform, data=None)

    def _process_quat_assembly(self, idx, data):
        """Create XtaQuatAssemblylform model instance from data.

        Incoming data format:
        assemblies:
        - chains: [A, B]
          assemblies_xtalforms: [0, 1, 2, 2, 3]

        No references to other models.
        """
        # TODO: unclear if there can be multiple or just one. the way the
        # data is structured it seems that there's only one.
        logger.debug("Creating QuatAssembly object: %d", idx)
        quat_assembly = QuatAssembly(
            chains=data["chains"],
            name="",  # TODO: name missing in yaml
        )
        try:
            quat_assembly.save()
        except IntegrityError as exc:
            msg = f"Failed to save QuatAssembly: {idx}"
            logger.error(msg)
            raise IntegrityError(msg) from exc

        return MetadataObjects(
            instance=quat_assembly,
            data=data["assemblies_xtalforms"],
        )

    def _process_canon_site_conf(self, canon_site, idx, data):
        """Create Xtalform model instance from data.

        Incoming data format:
        <idx>:
          name: <name>
          lig_ref: <lig_ref>
          residues:
            chain: <array[char[1]]>
            residue: <array[int]>
          members:
            dtag: <array[str], (protein codes)>
            chain: <array[char[1]]>
            residue: <array[int]>

        Unable to add references to:
        - SiteObservation (ref_site_observation)
        """
        logger.debug("Creating CanonSiteConf object: %s", data["name"])
        canon_site_conf = CanonSiteConf(
            name=data["name"],
            residues=data["residues"],
            canon_site=canon_site,
        )

        try:
            canon_site_conf.save()
        except IntegrityError as exc:
            msg = f"Failed to save CanonSiteConf: {data['name']}"
            logger.error(msg)
            raise IntegrityError(msg) from exc

        return MetadataObjects(
            instance=canon_site_conf,
            data=data["lig_ref"],
        )

    def _process_site_observation(
        self,
        experiment,
        compound,
        xtalform_site,
        canon_site_conf,
        chain,
        ligand,
        idx,
        data,
    ):
        """Create SiteObservation model instance from data.

        Incoming data format:
        <idx (apparently canon_site_conf id)>: {
          structure: <file path>,
          artefacts:  <file path>,
          event_map:  <file path>,
          x_map:  <file path>,
          pdb_apo:  <file path>,
          pdb_apo_solv:  <file path>,
          pdb_apo_desolv:  <file path>,
          ligand_mol:  <file path>,
          ligand_pdb:  <file path>,
          ligand_smiles: <smiles>,
        }
        """
        code = f"{experiment.code}_{chain}_{str(ligand)}_{str(idx)}"
        logger.debug("Creating SiteObservation object: %s", code)

        files = self._check_file_struct(self.target_root, data)

        mol_data = None
        with open(files["ligand_mol"], "r", encoding="utf-8") as f:
            mol_data = f.read()

        site_observation = SiteObservation(
            # Code for this protein (e.g. Mpro_Nterm-x0029_A_501_0)
            code=code,
            experiment=experiment,
            cmpd=compound,
            xtalform_site=xtalform_site,
            canon_site_conf=canon_site_conf,
            bound_file=str(files["structure"]),
            apo_solv_file=str(files["pdb_apo_solv"]),
            apo_desolv_file=str(files["pdb_apo_desolv"]),
            apo_file=str(files["pdb_apo"]),
            xmap_2fofc_file=str(files["x_map"]),
            event_file=str(files["event_map"]),
            artefacts_file=str(files["artefacts"]),
            trans_matrix_info="currently missing",
            pdb_header_file="currently missing",
            smiles=data["ligand_smiles"],
            seq_id=ligand,
            chain_id=chain,
            ligand_mol_file=mol_data,
        )

        try:
            site_observation.save()
        except IntegrityError as exc:
            msg = f"Failed to save SiteObservation: {site_observation.code}"
            logger.error(msg)
            raise IntegrityError(msg) from exc

        return MetadataObjects(
            instance=site_observation,
            data=None,
        )

    def _process_canon_site(self, idx, data):
        """Create CanonSite model instance from data.

        Incoming data format:
        <idx>:
          canon_site_ref_site: <int>
          canon_site_conf_sites: <array[int]>
          site_residues:
            chain: <array[char[1]]>
            residue:  <array[int]>
          site_members:
            dtag: <array[str], (protein codes)>
            chain: <array[char[1]]>
            residue: <array[int]>

        Unable to add references to:
        - QuatAssembly
        - CanonSiteConf (ref_conf_site)

        Need to add these later when they become available. (While
        QuatAssembly technically is, there's no straightforward reference
        to it.)

        """
        logger.debug("Creating CanonSite object: %d", idx)
        canon_site = CanonSite(
            name="",
            # TODO: unclear how to save residue info, open comment in spreadsheet
            residues=data["site_residues"],
            canon_site_id=idx,
        )
        try:
            canon_site.save()
        except IntegrityError as exc:
            msg = f"Failed to save CanonSite: {idx}"
            logger.error(msg)
            raise IntegrityError(msg) from exc

        return MetadataObjects(
            instance=canon_site,
            data=data,
        )

    def _process_xtalform_site(self, xtalform, canon_site, idx, data):
        """Create Xtalform model instance from data.

        Incoming data format:
        <idx>:
          xtalform_id: <int>
          canon_site_id: <int>
          lig_chain: A
          xtalform_members:
            dtag:  <array[str], (protein codes)>
            chain: <array[char[1]]>
            residue: <array[int]>

        Saves references to all other tables (Xtalform and CanonSite).
        """
        logger.debug("Creating XtalformSite object: %d", idx)
        xtalform_site = XtalformSite(
            xtalform=xtalform,
            canon_site=canon_site,
            lig_chain=data["lig_chain"],
            residues=data["xtalform_members"]["residue"],
            xtalform_site_id=idx,
        )

        try:
            xtalform_site.save()
        except IntegrityError as exc:
            msg = f"Failed to save Xtalform: {idx}"
            logger.error(msg)
            raise IntegrityError(msg) from exc

        return MetadataObjects(
            instance=xtalform_site,
            data=data["xtalform_members"],
        )

    def process_metadata(
        self, proposal_ref=None, contact_email=None, user_id=None, task=None
    ):
        """Extract model instances from metadata file and save them to db.

        If called from task, takes task as a parameter for status updates.

        The only method to be called outside of the class.
        """
        # TODO: this method is quite long and should perhaps be broken
        # apart. then again, it just keeps calling other methods to
        # create model instances, so logically it's just doing the
        # same thing.
        blocks = [
            "crystals",
            "assemblies",
            "xtalforms",
            "canon_sites",
            "conformer_sites",
            "xtalform_sites",
        ]
        try:
            (  # pylint: disable=unbalanced-tuple-unpacking
                crystals,
                assemblies,
                xtalforms,
                canon_sites,
                conformer_sites,
                xtalform_sites,
            ) = self._load_yaml_blocks(
                Path(self.target_root).joinpath(METADATA_FILE), blocks
            )
        except FileNotFoundError as exc:
            raise FileNotFoundError(exc.args[0]) from exc
        except KeyError as exc:
            raise KeyError(exc.args[0]) from exc

        try:
            task.update_state(
                state="PROCESSING",
                meta={
                    "description": "Processing metadata",
                },
            )
        except AttributeError:
            # no task passed to method, nothing to do
            pass

        target = get_create_target(self.target_name)

        # this is copied from get_create_projects. i'm not entirely
        # sure how it's used.
        visit = proposal_ref.split()[0]
        project = Project.objects.get_or_create(title=visit)[0]

        try:
            committer = get_user_model().objects.get(pk=user_id)
        except get_user_model().DoesNotExist:
            # TODO: need something here for testing
            committer = get_user_model().objects.get(pk=1)

        experiment_upload = ExperimentUpload(
            project=project,
            target=target,
            committer=committer,
            commit_datetime=timezone.now(),
            file=self.data_bundle,
        )
        experiment_upload.save()

        # saved in full, referenced in canonsite and xtalform
        # data has xtalform ids
        quat_assembly_objects = {}
        # unlike other blocks, this one is returned as list
        for idx, obj in enumerate(assemblies):
            try:
                quat_assembly_objects[idx] = self._process_quat_assembly(idx, obj)
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc

        experiment_objects = {}
        compound_objects = {}
        for prot_name, prot_data in crystals.items():
            # TODO: unclear if I need to save experiment that doesn't have
            # aligned crystal files section or should I move this beyond
            # continue as well.
            try:
                experiment_objects[prot_name] = self._process_experiment(
                    experiment_upload, prot_name, prot_data
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc

            cmpd_data = prot_data["crystallographic_files"]
            try:
                compound_objects[prot_name] = self._process_compound(
                    prot_name, cmpd_data
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc
            except KeyError:
                # this particular block doesn't have compound info
                # continue with the loop, nothing to do
                continue

        # save components manytomany to experiment
        # TODO: is it 1:1 relationship? looking at the meta_align it seems to be,
        # but why the m2m then?
        for comp_code, comp_meta in compound_objects.items():
            experiment = experiment_objects[comp_code].instance
            experiment.compounds.add(comp_meta.instance)

        # reindex quat assembly by xtalform
        quat_assembly_objects_by_xtalform = {
            xtalform: obj.instance
            for obj in quat_assembly_objects.values()
            for xtalform in obj.data
        }

        xtalform_objects = {}
        for idx, obj in xtalforms.items():
            try:
                xtalform_objects[idx] = self._process_xtalform(
                    experiment_objects[obj["xtalform_ref"]].instance,
                    quat_assembly_objects_by_xtalform[idx],
                    idx,
                    obj,
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc

        # up until this point, all objects are fully saved
        # missing fk's begin now

        canon_site_objects = {}
        for idx, obj in canon_sites.items():
            try:
                canon_site_objects[idx] = self._process_canon_site(idx, obj)
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc
        # NB! missing fk's:
        # - ref_conf_site
        # - quat_assembly

        # reindex canon sites by canon_sites_conf_sites
        # NB! this is also used below for ref_conf_site in canon_site
        canon_sites_by_conf_sites = {
            conf: obj.instance
            for obj in canon_site_objects.values()
            for conf in obj.data["canon_site_conf_sites"]
        }

        canon_site_conf_objects = {}
        for idx, obj in conformer_sites.items():
            try:
                canon_site_conf_objects[idx] = self._process_canon_site_conf(
                    canon_sites_by_conf_sites[idx],
                    idx,
                    obj,
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc
        # NB! missing fk's:
        # - site_observation

        xtalform_sites_objects = {}
        for idx, obj in xtalform_sites.items():
            try:
                xtalform_sites_objects[idx] = self._process_xtalform_site(
                    xtalform_objects[obj["xtalform_id"]].instance,
                    canon_site_objects[obj["canon_site_id"]].instance,
                    idx,
                    obj,
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc

        # now can update CanonSite with ref_conf_site and quat_assembly

        # xtalform_site references canon_site and xtalform, which means I
        # can now add ref to quat_assembly, via xtalform to canon_site
        canon_xtal_map = {
            k.instance.canon_site.pk: k.instance.xtalform
            for k in xtalform_sites_objects.values()
        }
        for val in canon_site_objects.values():
            val.instance.quat_assembly = canon_xtal_map[val.instance.pk].quat_assembly
            val.instance.ref_conf_site = canon_site_conf_objects[
                val.data["canon_site_ref_site"]
            ].instance
            val.instance.save()

        # canon site instances are now complete
        # still missing fk to site_observation in canon_site_conf

        # reindex xtalform site to grab for site observation
        xtalform_site_by_tag = {}
        for val in xtalform_sites_objects.values():
            for k in zip(val.data["dtag"], val.data["chain"], val.data["residue"]):
                xtalform_site_by_tag[k] = val.instance

        site_observation_objects = {}
        # TODO: would be nice to get rid of quadruple for
        for experiment_meta in experiment_objects.values():
            if experiment_meta.data is None:
                continue
            for chain, ligand in experiment_meta.data.items():
                for ligand, ligand_data in ligand.items():
                    for idx, obj in ligand_data.items():
                        try:
                            site_observation_objects[
                                (experiment_meta.instance.code, chain, ligand)
                            ] = self._process_site_observation(
                                experiment_meta.instance,
                                compound_objects[
                                    experiment_meta.instance.code
                                ].instance,
                                xtalform_site_by_tag[
                                    (experiment_meta.instance.code, chain, ligand)
                                ],
                                canon_site_conf_objects[idx].instance,
                                chain,
                                ligand,
                                idx,
                                obj,
                            )
                        except IntegrityError as exc:
                            raise IntegrityError(exc.args[0]) from exc

        # final remaining fk, attach reference site observation to canon_site_conf
        for val in canon_site_conf_objects.values():
            val.instance.ref_site_observation = site_observation_objects[
                (val.data["dtag"], val.data["chain"], val.data["residue"])
            ].instance

    @staticmethod
    def _check_file_struct(target_root, file_struct):
        """Check if file exists and if sha256 hash matches (if given).

        file struct can come in 2 configurations:
        {file_key: {file: <file_path>, sha265: <hash> [smiles: <smiles>]}, ...}
        or simply
        {file_key: <file path>}
        Detect which one and take appropriate action.
        """
        result = {}
        for key, value in file_struct.items():
            if isinstance(value, dict):
                try:
                    filename = value["file"]
                    # I guess I don't care if it's not given?
                    file_hash = value.get("sha256", None)
                    result[key] = TargetLoader._check_file(
                        target_root.joinpath(filename), file_hash
                    )
                except KeyError:
                    logger.warning("%s file info missing in %s!", key, METADATA_FILE)
                    result[key] = None
            elif isinstance(value, str):
                result[key] = TargetLoader._check_file(target_root.joinpath(value))
            else:
                # this is probably the list of panddas event files, don't
                # need them here
                # although.. should i validate them here nevertheless?
                # i'd have to do this on copy otherwise..
                pass

        return result

    @staticmethod
    def _check_file(file_path: Path, file_hash=None):
        """Check if file exist and compare with hash.

        file_dict -  {file: <file path>, sha256: <hash>}
        return: validated file paths or None
        """
        # TODO: unclear about the proper validation procedure. it
        # should ideally throw exception, not return Path or None.
        if file_path.is_file():
            if file_hash:
                if file_hash == TargetLoader._calculate_sha256(file_path):
                    return file_path
            else:
                return file_path
        return None

    # borrowed from SO
    @staticmethod
    def _calculate_sha256(filepath):
        sha256_hash = hashlib.sha256()
        with open(filepath, "rb") as f:
            # Read the file in chunks of 4096 bytes
            for chunk in iter(lambda: f.read(4096), b""):
                sha256_hash.update(chunk)
        return sha256_hash.hexdigest()


def load_target(
    data_bundle,
    proposal_ref=None,
    contact_email=None,
    user_id=None,
    task=None,
):
    # TODO: do I need to sniff out correct archive format?
    with TemporaryDirectory(dir=settings.MEDIA_ROOT) as tempdir:
        # for dev, skip extraction and open existing
        # tempdir = settings.MEDIA_ROOT
        target_loader = TargetLoader(data_bundle, tempdir)
        try:
            with tarfile.open(target_loader.bundle_path, "r") as archive:
                archive.extractall(target_loader.raw_data)
                # pass
        except FileNotFoundError as exc:
            msg = f"{data_bundle} file does not exist!"
            logger.error(msg)
            raise FileNotFoundError(msg) from exc

        try:
            with transaction.atomic():
                target_loader.process_metadata(
                    proposal_ref=proposal_ref,
                    contact_email=contact_email,
                    user_id=user_id,
                    task=task,
                )
        except FileNotFoundError as exc:
            raise FileNotFoundError(exc.args[0]) from exc
        except IntegrityError as exc:
            raise IntegrityError(exc.args[0]) from exc

        # data in db, all good, move the files to permanent storage

        # NB! removing existing data. Don't think this should happen
        # in production, but probably useful for testing.
        # TODO: remove when done
        path = (
            Path(settings.MEDIA_ROOT)
            .joinpath(TARGET_LOADER_DATA)
            .joinpath(target_loader.bundle_name)
        )
        if path.is_dir():
            import shutil

            logger.warning("Removing exiting data at %s", path)
            shutil.rmtree(path)

        target_loader.raw_data.rename(
            Path(settings.MEDIA_ROOT)
            .joinpath(TARGET_LOADER_DATA)
            .joinpath(target_loader.bundle_name)
        )

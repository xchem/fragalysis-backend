from pathlib import Path
from tempfile import TemporaryDirectory
import tarfile


from django.conf import settings


# from django.db import IntegrityError
from django.test import TestCase

# from viewer import target_loader

from viewer.target_loader import TargetLoader
from viewer.target_loader import MetadataObjects
from viewer.target_loader import ASSEMBLIES_FILE


from viewer.models import QuatAssembly


test_mpro_v1 = Path(__file__).absolute().parent.joinpath("Mpro-v1-zero.tgz")
test_mpro_v2 = Path(__file__).absolute().parent.joinpath("Mpro-v2-zero.tgz")


class LoaderTests(TestCase):
    tempdir = None
    target_loader = None

    @classmethod
    def setUpTestData(cls):
        cls.tempdir = TemporaryDirectory(dir=settings.MEDIA_ROOT)

        # set up target loader object
        # I'm not calling .name in live code, how is it working there??
        cls.target_loader = TargetLoader(test_mpro_v1, "lb-test", cls.tempdir.name)

        with tarfile.open(cls.target_loader.bundle_path, "r") as archive:
            archive.extractall(cls.target_loader.raw_data)

        # because I know where it is
        upload_root = Path(cls.target_loader.target_root).joinpath("upload_1")

        cls.assemblies = (
            cls.target_loader._load_yaml(  # pylint: disable=protected-access
                upload_root.joinpath(ASSEMBLIES_FILE)
            )
        )  # pylint: disable=protected-access

    @classmethod
    def tearDownClass(cls):
        cls.tempdir.cleanup()
        super().tearDownClass()

    def test__process_quat_assembly_positive(self):
        idx = next(iter(self.assemblies))
        data = self.assemblies[idx]

        result = self.target_loader._process_quat_assembly(  # pylint: disable=protected-access
            QuatAssembly.objects.none(), idx, data
        )  # pylint: disable=protected-access

        self.assertIsInstance(
            result, MetadataObjects, "Returned object is not of MetadataObjects type."
        )
        self.assertIsInstance(
            result.instance,
            QuatAssembly,
            "Returned instance is not of QuatAssembly type.",
        )

        # result.instance.pk = None
        # self.assertEqual(result.instance, quat_assembly, "Returned instance does not match the mock instance.")

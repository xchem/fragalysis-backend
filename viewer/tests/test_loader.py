from pathlib import Path


# from django.db import IntegrityError
from django.test import TestCase
# from viewer import target_loader

# from viewer.target_loader import TargetLoader
# from viewer.target_loader import MetadataObjects



# from viewer.models import QuatAssembly



assemblies_path = Path(__file__).absolute().parent.joinpath("assemblies.yaml")


class LoaderTests(TestCase):

    test_data = {'reference': '5rgs', 'biomol': 'A,B', 'chains': 'A, A(-x,y,-z)'}

    # def test__process_quat_assembly_positive(self):

    #     quat_assembly = QuatAssembly(
    #         chains=self.test_data["chains"],
    #         name="",
    #     )
    #     quat_assembly.save()

    #     idx = 1
    #     target_loader = TargetLoader()
    #     result = target_loader._process_quat_assembly(idx, self.test_data)

    #     print('quat_ass', repr(quat_assembly))
    #     print('resutl', result)


    #     self.assertIsInstance(result, MetadataObjects, "Returned object is not of MetadataObjects type.")
    #     self.assertEqual(result.instance, quat_assembly, "Returned instance does not match the mock instance.")
    #     # self.assertIsNone(result.data, "Returned data should be None.")

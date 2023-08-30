# from pathlib import Path

# from unittest.mock import MagicMock, patch

# from django.db import IntegrityError
# from django.test import TestCase

# from viewer.target_loader import TargetLoader
# from viewer.target_loader import MetadataObjects



# from viewer.models import QuatAssembly



# assemblies_path = Path(__file__).absolute().parent.joinpath("assemblies.yaml")


# class LoaderTests(TestCase):

#     test_data = {'reference': '5rgs', 'biomol': 'A,B', 'chains': 'A, A(-x,y,-z)'}

#     def test__process_quat_assembly_positive(self):

#         quat_assembly = QuatAssembly(
#             chains=self.test_data["chains"],
#             name="",
#         )
#         quat_assembly.save()

#         idx = 1
#         with patch("viewer.target_loader.TargetLoader") as MockTargetLoader:
#             instance_mock = MagicMock()
#             # MockTargetLoader.return_value = instance_mock
#             # instance_mock.save.return_value = None
#             result = instance_mock._process_quat_assembly(idx, self.test_data)

#         print('instance mock', repr(instance_mock))
#         print('mockquat', MockTargetLoader)
#         print('quat_ass', repr(quat_assembly))
#         print('resutl', result)


#         self.assertIsInstance(result, MetadataObjects, "Returned object is not of MetadataObjects type.")
#         self.assertEqual(result.instance, instance_mock, "Returned instance does not match the mock instance.")
#         self.assertIsNone(result.data, "Returned data should be None.")


# from django.db import IntegrityError
# from django.test import TestCase
# from unittest.mock import MagicMock, patch
# from your_module import create_quat_assembly, QuatAssembly, MetadataObjects

# class YourTestClass(TestCase):
#     def test_create_quat_assembly_positive(self):
#         test_data = {
#             "reference": "TestAssembly",
#             "biomol": "biomol_data",
#             "chains": "chain_info"
#         }

#         with patch("your_module.QuatAssembly") as MockQuatAssembly:
#             instance_mock = MagicMock()
#             MockQuatAssembly.return_value = instance_mock
#             instance_mock.save.return_value = None

#             result = create_quat_assembly(test_data)

#         self.assertIsInstance(result, MetadataObjects, "Returned object is not of MetadataObjects type.")
#         self.assertEqual(result.instance, instance_mock, "Returned instance does not match the mock instance.")
#         self.assertIsNone(result.data, "Returned data should be None.")

#         MockQuatAssembly.assert_called_once_with(
#             chains=test_data["chains"],
#             name="",
#         )
#         instance_mock.save.assert_called_once()

#     def test_create_quat_assembly_negative(self):
#         test_data = {
#             "reference": "TestAssembly",
#             "biomol": "biomol_data",
#             "chains": "chain_info"
#         }

#         with patch("your_module.QuatAssembly") as MockQuatAssembly:
#             instance_mock = MagicMock()
#             MockQuatAssembly.return_value = instance_mock
#             instance_mock.save.side_effect = IntegrityError()

#             with self.assertRaises(IntegrityError):
#                 create_quat_assembly(test_data)

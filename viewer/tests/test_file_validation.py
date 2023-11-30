from pathlib import Path
import copy

from unittest import TestCase

from viewer.target_loader import TargetLoader


test_file1_path = Path(__file__).absolute().parent.joinpath("hash_test_file1.txt")
test_file1_hash = "35360b39bff52650367f7d2c08e29fa00a63a374c777b4055c75e4f97b173271"

test_file2_path = Path(__file__).absolute().parent.joinpath("hash_test_file2.txt")
test_file2_hash = "f2b6d99b58c4b32966779533b0987a769c11f0cfd9be4f20520a97032e41cb23"


file_struct_flat = {
    "file1": str(test_file1_path),
    "file2": str(test_file2_path),
}

file_struct_nested = {
    "file1": {
        "file": str(test_file1_path),
        "sha256": test_file1_hash,
    },
    "file2": {
        "file": str(test_file2_path),
        "sha256": test_file2_hash,
    },
}


class FileValidationTests(TestCase):
    def test__calculate_sha256_positive(self):
        calculated_hash = (
            TargetLoader._calculate_sha256(  # pylint: disable=protected-access
                test_file1_path
            )
        )  # pylint: disable=protected-access

        self.assertEqual(
            calculated_hash,
            test_file1_hash,
            "Hashes do not match for the positive test case.",
        )

    def test__calculate_sha256_negative(self):
        incorrect_hash = "imagine if this were the actual hash"
        calculated_hash = (
            TargetLoader._calculate_sha256(  # pylint: disable=protected-access
                test_file1_path
            )
        )  # pylint: disable=protected-access

        self.assertNotEqual(
            calculated_hash,
            incorrect_hash,
            "Hashes should not match for the negative test case.",
        )

    def test__check_file_existence(self):
        file_path = Path(test_file1_path)
        result = TargetLoader._check_file(file_path)  # pylint: disable=protected-access
        self.assertTrue(
            result, "File existence check failed for the positive test case."
        )

    def test__check_nonexistent_file(self):
        invalid_file_path = Path("path/to/nonexistent/file.txt")
        result = TargetLoader._check_file(  # pylint: disable=protected-access
            invalid_file_path
        )  # pylint: disable=protected-access
        self.assertFalse(
            result, "File existence check succeeded for the nonexistent file case."
        )

    def test__check_file_struct_flat_positive(self):
        result = TargetLoader._check_file_struct(  # pylint: disable=protected-access
            Path(test_file1_path.root), file_struct_flat
        )

        self.assertEqual(
            result,
            file_struct_flat,
            "File structure check failed for the positive test case.",
        )

    def test__check_file_struct_flat_incomplete_positive(self):
        file_struct = copy.deepcopy(file_struct_flat)
        del file_struct["file2"]

        result = TargetLoader._check_file_struct(  # pylint: disable=protected-access
            Path(test_file1_path.root), file_struct
        )

        expected_result = {
            "file1": str(test_file1_path),
        }

        self.assertEqual(
            result,
            expected_result,
            "File structure check failed for the positive test case.",
        )

    def test_check_file_struct_nested_positive(self):
        result = TargetLoader._check_file_struct(  # pylint: disable=protected-access
            Path(test_file1_path.root), file_struct_nested
        )

        self.assertEqual(
            result,
            file_struct_flat,
            "File structure check failed for the positive test case.",
        )

    def test_check_file_struct_nested_incomplete_positive(self):
        file_struct = copy.deepcopy(file_struct_nested)
        file_struct["file2"]["sha256"] = "incorrect hash"

        expected_result = copy.deepcopy(file_struct_flat)
        del expected_result["file2"]

        result = TargetLoader._check_file_struct(  # pylint: disable=protected-access
            Path(test_file1_path.root), file_struct
        )

        self.assertEqual(
            result,
            expected_result,
            "File structure check failed for the positive test case.",
        )

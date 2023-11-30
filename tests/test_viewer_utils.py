import filecmp
import os
from unittest.mock import Mock, patch

from django.test import SimpleTestCase, tag

from viewer.utils import (
    create_squonk_job_request_url,
    create_media_sub_directory,
    delete_media_sub_directory,
    add_prop_to_sdf,
    clean_filename,
    get_https_host,
)


@tag("nodb")
class ViewerUtilsTestCase(SimpleTestCase):
    def test_create_squonk_job_request_url(self):
        result = create_squonk_job_request_url("instance-0000")
        self.assertEqual(result, "data-manager-ui/results/instance/instance-0000")

    def test_create_media_sub_directory(self):
        result = create_media_sub_directory("blob")
        self.assertEqual(result, "/code/media/blob")
        test_rel_path = result[6:]
        self.assertTrue(os.path.exists(test_rel_path))
        os.rmdir(test_rel_path)

    def test_delete_media_sub_directory(self):
        result = create_media_sub_directory("blobby")
        self.assertEqual(result, "/code/media/blobby")
        test_rel_path = result[6:]
        self.assertTrue(os.path.exists(test_rel_path))
        delete_media_sub_directory("blobby")
        test_rel_path = result[6:]
        self.assertFalse(os.path.exists(test_rel_path))

    def test_add_prop_to_sdf(self):
        expected_file = "./tests/_test.sdf"
        self.assertFalse(os.path.isfile(expected_file))
        add_prop_to_sdf(
            "tests/test_data/viewer-utils-test-in.sdf",
            expected_file,
            {"TransFSScore": "0.115601"},
        )
        self.assertTrue(os.path.isfile(expected_file))
        self.assertTrue(
            filecmp.cmp("tests/test_data/viewer-utils-test-out.sdf", expected_file)
        )
        os.remove("./tests/_test.sdf")

    def test_clean_filename_a(self):
        result = clean_filename("./media/sdfs/Mpro-x3351_0A_rtEVbqf.sdf")
        self.assertEqual(result, "Mpro-x3351_0A.sdf")

    def test_clean_filename_b(self):
        result = clean_filename("./media/sdfs/Mpro-x3351_0A.sdf")
        self.assertEqual(result, "Mpro-x3351_0A.sdf")

    def test_clean_filename_c(self):
        result = clean_filename("Mpro-x3351_0A.sdf")
        self.assertEqual(result, "Mpro-x3351_0A.sdf")

    def test_get_https_host(self):
        request_mock = Mock()
        request_mock.get_host.return_value = "example.com"
        result = get_https_host(request_mock)
        self.assertEqual(result, "https://example.com")

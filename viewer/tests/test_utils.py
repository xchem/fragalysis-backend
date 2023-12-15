import filecmp
import os
from unittest.mock import Mock

from django.test import SimpleTestCase

from viewer.utils import (
    add_prop_to_sdf,
    clean_filename,
    create_media_sub_directory,
    create_squonk_job_request_url,
    delete_media_sub_directory,
    get_https_host,
    is_url,
    word_count,
)


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

    def test_is_url(self):
        self.assertTrue(is_url("https://example.com"))
        self.assertTrue(is_url("http://example.com"))
        self.assertTrue(is_url("ftp://example.com"))
        self.assertTrue(is_url("ftps://example.com"))
        self.assertTrue(is_url("sftp://example.com"))
        self.assertTrue(is_url("ssh://example.com"))
        self.assertTrue(is_url("file://example.com"))
        self.assertTrue(is_url("mailto://example.com"))
        self.assertTrue(is_url("telnet://example.com"))
        self.assertTrue(is_url("gopher://example.com"))
        self.assertTrue(is_url("wais://example.com"))
        self.assertTrue(is_url("ldap://example.com"))
        self.assertTrue(is_url("news://example.com"))
        self.assertTrue(is_url("nntp://example.com"))
        self.assertTrue(is_url("prospero://example.com"))
        # SOme failures...
        self.assertFalse(is_url("/data/blob.html"))
        self.assertFalse(is_url(532))  # type: ignore[arg-type]
        self.assertFalse(is_url(None))

    def test_word_count(self):
        self.assertEquals(2, word_count("Hello world"))
        self.assertEquals(1, word_count("Hello"))
        self.assertEquals(0, word_count(""))
        self.assertEquals(0, word_count(None))

"""A Django storage class for S3, built on boto3."""
from .storage import S3Storage, S3UploadError

__all__ = ["S3Storage", "S3UploadError"]

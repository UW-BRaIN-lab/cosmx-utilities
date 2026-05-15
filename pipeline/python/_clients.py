"""Shared boto3 client factories for the pipeline.

Two storage backends:
  - AWS S3 (source of truth): credentials resolved from AWS_SOURCE_*, or
    AWS_SOURCE_PROFILE, or boto3's default credential chain.
  - Kopah (working storage):  credentials resolved from KOPAH_ACCESS_KEY_ID /
    KOPAH_SECRET_ACCESS_KEY against KOPAH_ENDPOINT_URL.

Keeping these in one module prevents the AWS_* env vars set by the klone
`kopah_brain` shell helper from leaking into the AWS source client.
"""

from __future__ import annotations

import os

import boto3
from botocore.client import Config


def make_source_client():
    """AWS S3 client, with credentials resolved in priority order:

    1. AWS_SOURCE_ACCESS_KEY_ID + AWS_SOURCE_SECRET_ACCESS_KEY (explicit)
    2. AWS_SOURCE_PROFILE (named profile from ~/.aws/credentials)
    3. boto3's default credential chain
    """
    access_key = os.environ.get("AWS_SOURCE_ACCESS_KEY_ID")
    secret_key = os.environ.get("AWS_SOURCE_SECRET_ACCESS_KEY")
    if access_key and secret_key:
        return boto3.client(
            "s3",
            aws_access_key_id=access_key,
            aws_secret_access_key=secret_key,
        )
    profile = os.environ.get("AWS_SOURCE_PROFILE")
    if profile:
        return boto3.Session(profile_name=profile).client("s3")
    return boto3.client("s3")


def make_kopah_client(endpoint_url: str, access_key: str, secret_key: str):
    # request_checksum_calculation / response_checksum_validation set to
    # "when_required" disable boto3 1.36+'s default-on x-amz-checksum-* headers
    # for uploads. Ceph RGW (which Kopah runs on) rejects those headers during
    # multipart UploadPart with XAmzContentSHA256Mismatch.
    return boto3.client(
        "s3",
        endpoint_url=endpoint_url,
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
        config=Config(
            signature_version="s3v4",
            s3={"addressing_style": "path"},
            request_checksum_calculation="when_required",
            response_checksum_validation="when_required",
        ),
    )

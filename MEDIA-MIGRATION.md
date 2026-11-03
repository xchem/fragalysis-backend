# Moving a stack's media to S3 — runbook

This runbook moves an existing stack's media files from the filesystem (NFS,
`MEDIA_STORAGE=FILESYSTEM`) to an S3 bucket (`MEDIA_STORAGE=S3`), for example STFC
Echo. It is for the people who deploy and operate a stack. The design is in
issue #1038, and `ARCHITECTURE.md` (Media storage) describes the code.

The layout does not change: an object's key is the file's path relative to the media
directory (under `MEDIA_S3_LOCATION`, if it is set). So a rollback is a copy in the
other direction.

## Before you start

1. **The bucket.** Make the bucket on the S3 service, and get credentials that can
   read, write, list and delete in it. Keep the credentials in a Kubernetes Secret.
2. **The settings.** Add these to the stack's deployment (Ansible/Helm), for the
   stack and the Celery workers. Keep `MEDIA_STORAGE=FILESYSTEM` for now.

   | Variable | Value |
   |---|---|
   | `MEDIA_S3_BUCKET_NAME` | The bucket name |
   | `MEDIA_S3_ENDPOINT_URL` | i.e. `https://s3.echo.stfc.ac.uk`. In Kubernetes, the host must be a fully qualified domain name: the nginx resolver does not use search domains. |
   | `MEDIA_S3_ACCESS_KEY_ID`, `MEDIA_S3_SECRET_ACCESS_KEY` | From the Secret |
   | `MEDIA_S3_REGION_NAME` | If the service needs one |
   | `MEDIA_S3_LOCATION` | Optional: a key prefix, i.e. `media` |
   | `MEDIA_TRANSFER_THREADS` | Optional (default `16`) |
   | `MEDIA_TMP_MAX_AGE_H` | Optional (default `48`) |
   | `NGINX_RESOLVER` | Optional: only if `/etc/resolv.conf` does not give a usable DNS server |

3. **Scratch space.** In S3 mode, loads, downloads and computed-set uploads work in
   Pod-local scratch (`/tmp`). Give the stack and worker Pods an emptyDir for `/tmp`
   of about twice the largest bundle or download.
4. **Network.** The stack and worker Pods must reach the S3 endpoint (HTTPS, 443).
5. **Check the S3 service.** The code uses `CopyObject` (to move a file), batch
   `DeleteObjects`, presigned GET URLs, `Range` requests and multipart uploads. On
   Echo, all except `CopyObject` are confirmed (see #1038). Confirm `CopyObject` with
   the bucket's credentials:

   ```bash
   aws --endpoint-url https://s3.echo.stfc.ac.uk s3 cp ./test.txt s3://<bucket>/test.txt
   aws --endpoint-url https://s3.echo.stfc.ac.uk s3 cp s3://<bucket>/test.txt s3://<bucket>/test-copy.txt
   aws --endpoint-url https://s3.echo.stfc.ac.uk s3 rm s3://<bucket>/test.txt
   aws --endpoint-url https://s3.echo.stfc.ac.uk s3 rm s3://<bucket>/test-copy.txt
   ```

## The steps

Run `copy_media` in a backend container (i.e. in the stack Pod) that has the media
directory mounted and the `MEDIA_S3_*` settings. It works whatever `MEDIA_STORAGE`
is.

1. **Release the code in FILESYSTEM mode.** Migration `0174` changes
   `ComputedSet.written_sdf_filename` to storage names. Check that the stack works as
   before.
2. **Dry run.** This copies nothing. It gives the files and bytes for each top-level
   prefix, which are also the numbers for the bucket's quota:

   ```bash
   python manage.py copy_media --dry-run
   ```

3. **Copy while the stack runs.** This can take hours. If it stops, run it again: it
   copies only the files that are missing (or have a different size) in the bucket.

   ```bash
   python manage.py copy_media
   ```

4. **The maintenance window.**
   1. Stop user access to the stack (i.e. a maintenance page at the ingress).
   2. Scale the Celery workers to 0, so that no upload or task changes the media.
   3. Copy the changes since step 3:

      ```bash
      python manage.py copy_media
      ```

   4. Verify the copy. It must pass:

      ```bash
      python manage.py copy_media --verify
      ```

   5. Check the database references. It must pass. Read its list of files that are
      *missing on the source too*: they were missing on NFS before the migration.

      ```bash
      python manage.py copy_media --check-references
      ```

   6. Set `MEDIA_STORAGE=S3` for the stack and the workers, and restart them. The
      start-up check (`viewer.E003`) stops the stack if it cannot use the bucket.
   7. Test the stack:
      - Open a target: the 3D viewer loads its PDBs and maps.
      - In the browser's developer tools, the media requests go to `/media/...` on the
        stack's own host, and no response has an `X-Media-Url` header.
      - Download structures: the size is shown (a HEAD request), then the archive
        downloads.
      - Upload a computed set and a (test) target, and delete the test target.
   8. Open user access.

5. **Keep NFS** mounted read-only for the rollback period. Then remove the media
   PersistentVolumeClaim.

What `copy_media` does not copy:

- `tmp/` (use `--include-tmp` to copy it);
- `squonk_upload/` (the Squonk2 integration is deprecated);
- the media deletion watcher's files;
- download archives other than those of static links and of dynamic links that have
  not expired. The others cannot be fetched.

It reports, and does not copy, symlinks, names that are not UTF-8, and names whose
S3 key would be longer than 1024 bytes.

## Rollback

- **Back to the filesystem, on the same release.** Stop user access and the workers.
  Make the NFS volume writable. Copy the files that were made after the cutover back
  to it, then set `MEDIA_STORAGE=FILESYSTEM` and restart:

  ```bash
  python manage.py copy_media --direction s3-to-fs
  ```

- **Back to a release before the S3 work.** Do the step above, then reverse migration
  `0174` before you deploy the older release:

  ```bash
  python manage.py migrate viewer 0173
  ```

Files that were deleted in S3 after the cutover are still on NFS. They are not
referenced, and `cleanup_media` can remove them.

## After the migration

- Remove the media PersistentVolumeClaim (step 5).
- Agree with the S3 service's operators (i.e. STFC) on backups or versioning for the
  bucket, and on access logs.
- Watch the bucket's usage against its quota.

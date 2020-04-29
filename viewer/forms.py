from django import forms
import zipfile
from cStringIO import StringIO

CHOICES = [
    (0, 'validate'),
    (1, 'upload'),
]

class CSetForm(forms.Form):
    target_name = forms.CharField(label='Target', max_length=100)
    sdf_file = forms.FileField(label='All compounds sdf (.sdf)')
    pdb_zip = forms.FileField(required=False, label='PDB files (.zip)')
    submit_choice = forms.CharField(widget=forms.RadioSelect(choices=CHOICES))
    upload_key = forms.CharField(label='Upload Key')

    # def clean_zipfile(self):
    #     if self.zip_file.get('content-type') != 'application/zip':
    #         msg = 'File upload must be a valid ZIP archive.'
    #         raise forms.ValidationError(msg)
    #     else:
    #         try:
    #             zip = zipfile.ZipFile(StringIO(self.zip_file['content']))
    #         except:
    #             raise forms.ValidationError("Could not unzip file.")
    #         bad_file = zip.testzip()
    #         zip.close()
    #         del zip
    #         if bad_file:
    #             raise forms.ValidationError(msg)
    #     return self.zip_file  # Return the clean zip_file

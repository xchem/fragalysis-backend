from django import forms
import zipfile
from io import StringIO

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


class UploadKeyForm(forms.Form):
    contact_email = forms.EmailField(widget=forms.TextInput(attrs={'class':'form-control', 'autocomplete':'off'}), required=True)

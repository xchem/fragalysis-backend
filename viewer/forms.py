from django import forms

class CSetForm(forms.Form):
    target_name = forms.CharField(label='Target', max_length=100)
    sdf_file = forms.FileField()
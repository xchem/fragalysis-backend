from django import forms

CHOICES = [
    (0, 'validate'),
    (1, 'upload'),
]

class CSetForm(forms.Form):
    target_name = forms.CharField(label='Target', max_length=100)
    sdf_file = forms.FileField()
    submit_choice = forms.CharField(widget=forms.RadioSelect(choices=CHOICES))
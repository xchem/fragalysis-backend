from django import forms

class CSetForm(forms.Form):
    target_name = forms.CharField(label='Target', max_length=100)
    sdf_file = forms.FileField()
    submit_choice = forms.ChoiceField(widget=forms.RadioSelect, choices=[(0, 'validate'), (1, 'upload')])
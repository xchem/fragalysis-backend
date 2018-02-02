from .models import Protein
from django import forms

class PDBForm(forms.ModelForm):
    class Meta:
        model = Protein
        fields = ('file', )
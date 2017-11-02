from django.http import HttpResponse
import json
from django.shortcuts import render

def display(request):
    return render(request, 'viewer/display.html', {})
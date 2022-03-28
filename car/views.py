import os
from django.views import View
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from rest_framework.decorators import action
from django.http import JsonResponse
from django.http import JsonResponse
from django.conf import settings
from celery.result import AsyncResult

from viewer.tasks import check_services
from .tasks import (
    validateFileUpload,
    uploadManifoldReaction,
    uploadCustomReaction,
)
import pandas as pd


def save_tmp_file(myfile):
    name = myfile.name
    path = default_storage.save("tmp/" + name, ContentFile(myfile.read()))
    tmp_file = str(os.path.join(settings.MEDIA_ROOT, path))
    return tmp_file


class UploadProject(View):
    @action(methods=['post'], detail=False)
    def uploadproject(self, request, pk=None):
        check_services()
        print(request)
        project_info = {}
        project_info["submittername"] = request.data["submitter_name"]
        project_info["submitterorganisation"] = request.data["submitter_organisation"]
        project_info["submitteremail"] = request.data["submitter_email"]
        validate_choice = request.data["validate_choice"]
        API_choice = request.data["API_choice"]

        csvfile = request.FILES["csv_file"]
        tmp_file = save_tmp_file(csvfile)

        if str(validate_choice) == "0":

                if str(API_choice) == "1":
                    task = validateFileUpload.delay(
                        csv_fp=tmp_file, validate_type="custom-chem"
                    )
                    
                if str(API_choice) == "2":
                    task = validateFileUpload.delay(
                        csv_fp=tmp_file, validate_type="combi-custom-chem"
                    )

                else:
                    task = validateFileUpload.delay(
                        csv_fp=tmp_file, validate_type="retro-API"
                    )

        if str(validate_choice) == "1":
                if str(API_choice) == "0":
                    task = (
                        validateFileUpload.s(
                            csv_fp=tmp_file,
                            validate_type="retro-API",
                            project_info=project_info,
                            validate_only=False,
                        )
                        | uploadManifoldReaction.s()
                    ).apply_async()
            
                if str(API_choice) == "1":
                    task = (
                        validateFileUpload.s(
                            csv_fp=tmp_file,
                            validate_type="custom-chem",
                            project_info=project_info,
                            validate_only=False,
                        )
                        | uploadCustomReaction.s()
                    ).apply_async()

                if str(API_choice) == "2":
                    task = (
                        validateFileUpload.s(
                            csv_fp=tmp_file,
                            validate_type="combi-custom-chem",
                            project_info=project_info,
                            validate_only=False,
                        )
                        | uploadCustomReaction.s()
                    ).apply_async()

        data = {"task_id": task.id}
        return JsonResponse(data=data)
    
    @action(detail=False, methods=['get'])
    def gettaskstatus(self, request, pk=None):
        task_id = self.request.GET.get('task_id', None)
        if task_id:
            task = AsyncResult(task_id)            
            if task.status == "FAILURE":
                data = {"task_status": task.status, "traceback": str(task.traceback)}
                return JsonResponse(data)

            if task.status == "SUCCESS":
                results = task.get()
                validate_dict = results[0]
                validated = results[1]
            if validated:
                data = {"task_status": task.status, "task_summary": "Success"}

                return JsonResponse(data)

            if not validated:
                pd.set_option("display.max_colwidth", None)
                table = pd.DataFrame.from_dict(validate_dict)
                html_table = table.to_html()
                html_table += (
                    """<p> Your data was <b>not</b> validated. The table above shows errors</p>"""
                )

                data = {"task_status": task.status, "task_summary": html_table}

                return JsonResponse(data)
                
            if task.status == "PENDING":
                data = {"task_status": task.status}
                return JsonResponse(data)
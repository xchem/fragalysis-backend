from rest_framework.test import APIClient, APITestCase
from car.createmodels import createBatchModel, createProjectModel

from car.models import Project

class APICreateProjectTestCase(APITestCase):
    def setUp(self) -> None:
        self.project_id = createProjectModel(project_info={"projectname": "Test", "submittername": "Tester", "submitterorganisation": "Rest Framework", "proteintarget": "MID2"})
        
        # for batchtag, group in grouped_targets:
        #     self.batch_id = createBatchModel(project_id=self.project_id, batchtag=batchtag)
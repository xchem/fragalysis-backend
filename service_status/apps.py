from django.apps import AppConfig


class ServiceStatusConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'service_status'

    def ready(self):
        # dummy import needed because otherwise tasks aren't being registered
        import service_status.services  # pylint: disable=unused-import

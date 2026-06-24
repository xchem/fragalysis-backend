from django.apps import AppConfig


class ServiceStatusConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'service_status'

    def ready(self):
        # Import ensures service functions have their _is_service_query marker set
        import service_status.services  # pylint: disable=unused-import

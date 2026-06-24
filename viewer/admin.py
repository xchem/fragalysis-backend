from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from django.contrib.auth.models import User

from viewer.models import UserRole


class UserRoleInline(admin.TabularInline):
    model = UserRole.users.through
    extra = 1
    verbose_name = "role"
    verbose_name_plural = "roles"


@admin.register(UserRole)
class UserRoleAdmin(admin.ModelAdmin):
    list_display = ("name",)
    search_fields = ("name",)
    ordering = ("name",)


class FragalysisUserAdmin(UserAdmin):
    inlines = (UserRoleInline,)


admin.site.unregister(User)
admin.site.register(User, FragalysisUserAdmin)

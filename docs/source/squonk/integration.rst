Squonk Integration
==================

The `Squonk`_ **Data Manager** and associated services, developed by `Informatics Matters`_,
provides a novel, easy to use, web based, data centric workflow environment in which
scientists can execute scientific workflows using open source and commercial tools from
multiple sources such as RDKit, Chemistry Development Kit (CDK), ChemAxon.

A RESTful API is also available for accessing the data and services **Squonk** provides.

Fragalysis has been adapted so that it can utilise **Squonk** services via its API
and Graphical User Interface.

Squonk Installation
-------------------

Rather than use the commercial **Squonk** service a *custom* **Squonk** installation
is typically provided, running in the same cluster as the Fragalysis Stack that
you wish to use. The installation will consist of the **Squonk** *Data Manager*
it accounting service, the *account server*. The *Data Manager* and *Account Server*
will be configured to use the keycloak authentication server also used by Fragalysis.

The administrator of these applications will need to ensure the following: -

Fragalysis Squonk User
^^^^^^^^^^^^^^^^^^^^^^

An administrative user is made available for use exclusively by Fragalysis. This is
the user Fragalysis will use to automate the orchestration of objects in **Squonk** -
the creation of *Units*, *Products*, *Projects* and Jobs (*Instances*). The user name
is typically `fragalysis`.

.. note::
    If more than one Fragalysis Stack is configured to use one **Squonk** installation
    the user can be shared between them - you don't necessarily need a separate
    `fragalysis` user for each Fragalysis Stack.

Fragalysis Squonk Organisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An *Organisation* needs to be created that will be the root of all objects created
by Fragalysis (its *Units*, *Products* and *Projects*). You don't need to know
how these objects relate to Fragalysis activity, only that an *Organisation* must
be provisioned.

.. note::
    The Organisation is Stack-specific - **EVERY** Fragalysis Stack must be assigned
    a unique organisation. The *Staging* and *Production* stacks currently share
    access to One **Squonk**, and they each do this using separate *Organisations*.

Squonk User Roles
^^^^^^^^^^^^^^^^^

Users that use Fragalysis, and also need access to **Squonk** services,
will need to be assigned the ``data-manager-user`` and ``account-server-user`` **Roles**
in the corresponding keycloak service.

Squonk Integration
------------------

With the **Squonk** installation in place, the Fragalysis Stack needs be configured
using the following Ansible playbook variables, which translate to environment variables
available to the Fragalysis Stack container: -

.. list-table:: Stack Playbook Variables
   :widths: 30 70
   :header-rows: 1
   
   * - Variable
     - Description
   * - ``stack_oidc_as_client_id``
     - The client ID of the Account Server in the keycloak authentication server.
   * - ``stack_oidc_dm_client_id``
     - The client ID of the Data Manager in the keycloak authentication server.
   * - ``stack_squonk2_ui_url``
     - The URL of the installed Squonk UI, typically ``https://squonk.example.com/data-manager-ui``
   * - ``stack_squonk2_dmapi_url``
     - The URL of the installed Squonk Data Manager API, typically ``https://squonk.example.com/data-manager-api``
   * - ``stack_squonk2_asapi_url``
     - The URL of the installed Squonk Account Server API, typically ``https://squonk.example.com/account-server-api``
   * - ``stack_squonk2_org_owner``
     - The name of the Fragalysis user that will be used to create objects in Squonk.
   * - ``stack_squonk2_org_owner_password``
     - The password of the Fragalysis user that will be used to create objects in Squonk.
   * - ``stack_squonk2_org_uuid``
     - The UUID of the Squonk Organisation under which *Units*, *Products* and *Projects*
       will be created. This must be unique for each Fragalysis Stack.
   * - ``stack_squonk2_product_flavour``
     - The *flavour* of the Squonk *Products* that will be created.
       This must be one of ``BRONZE``, ``SILVER`` or ``GOLD``.
   * - ``stack_squonk2_slug``
     - The *slug* used to create objects in Squonk *Product*. This is used to
       identify the Fragalysis Stack that created the object and must be
       unique for each Fragalysis Stack. It is a short string (up to 10 characters).
       The staging stack might use ``staging`` and the production stack might use
       ``production``.
   * - ``stack_squonk2_unit_billing_day``
     - The day of the month on which the *Units* will be billed. This must be
       an integer between 1 and 28.

.. warning::
    The ``stack_squonk2_org_uuid``, which must be unique for each Fragalysis Stack,
    cannot currently be changed once the stack has been launched. It is the
    responsibility of the administrator to ensure that the UUID is unique and the
    variable is correctly set prior to launching the stack.

Squonk Model
------------

The following Fragalysis ``viewer.models`` are used by the Fragalysis Stack to manage
access to and record the **Squonk** objects that are created: -

.. autoclass:: viewer.models.Squonk2Org
    :members:

.. autoclass:: viewer.models.Squonk2Unit
    :members:

.. autoclass:: viewer.models.Squonk2Project
    :members:

.. autoclass:: viewer.models.JobFileTransfer
    :members:

.. autoclass:: viewer.models.JobRequest
    :members:

Squonk Views
------------

.. autoclass:: viewer.views.JobConfigView
    :members:

.. autoclass:: viewer.views.JobFileTransferView
    :members:

.. autoclass:: viewer.views.JobRequestView
    :members:

.. autoclass:: viewer.views.JobCallBackView
    :members:

.. autoclass:: viewer.views.JobAccessView
    :members:

Squonk2Agent Class
------------------

The main interactions with **Squonk** are handled by the ``Squonk2Agent`` class
in the ``viewer`` package and the `Squonk2 Python Client`_, which provides API access
to **Squonk** through its ``DmApi`` and ``AsApi`` classes.

.. autoclass:: viewer.squonk2_agent.Squonk2Agent
    :members:

.. _Informatics Matters: https://www.informaticsmatters.com/
.. _Squonk: https://squonk.it/
.. _Squonk2 Python Client: https://github.com/InformaticsMatters/squonk2-python-client

"""
discourse.py Functions for to send and retrieve posts to/from the discourse platform.
Discourse is not considered available if the DISCOURSE_API_KEY is not set.
"""
import logging
import os
from typing import Tuple

import pydiscourse
from django.conf import settings

from viewer.models import DiscourseCategory, DiscourseTopic

logger = logging.getLogger(__name__)


def _get_client():
    """Create Discourse client for the configured fragalysis user"""
    assert settings.DISCOURSE_API_KEY
    return pydiscourse.DiscourseClient(
        settings.DISCOURSE_HOST,
        api_username=settings.DISCOURSE_USER,
        api_key=settings.DISCOURSE_API_KEY,
    )


def _get_user_client(user):
    """Create Discourse client for user.
    The code would normally have checked that the user exists before calling this.
    """
    assert settings.DISCOURSE_API_KEY
    return pydiscourse.DiscourseClient(
        settings.DISCOURSE_HOST,
        api_username=user.username,
        api_key=settings.DISCOURSE_API_KEY,
    )


def get_user(client, username) -> Tuple[bool, str, int]:
    """Call discourse API users to retrieve user by username."""
    assert client

    logger.info('+ discourse.get_user(%s)', username)

    error: bool = False
    error_message: str = ''
    user: int = 0

    try:
        user_record = client.user(username)
    except Exception as e:
        # User is not known in Discourse or there is a failure accessing Discourse.
        logger.error('discourse.get_user', exc_info=e)
        error = True
        if settings.DISCOURSE_HOST:
            error_message = (
                'Error validating user in Discourse. If this is your first post '
                + 'please log on to '
                + 'Discourse once to create a User. URL is: '
                + settings.DISCOURSE_HOST
            )
        else:
            error_message = 'Please check Discourse Host parameter for Fragalysis'
    else:
        user = user_record['id']

    logger.info(
        '- discourse.get_user(%s) error=%s error_message="%s"',
        username,
        error,
        error_message,
    )
    return error, error_message, user


def create_category(
    client,
    category_name,
    parent_name,
    category_colour='0088CC',
    category_text_colour='FFFFFF',
):
    """Call discourse API categories.json to create a new (sub)category."""
    assert client

    logger.info('+ discourse.create_category(%s)', category_name)

    category = client.create_category(
        category_name,
        color=category_colour,
        text_color=category_text_colour,
        parent=parent_name,
    )
    topic_url = os.path.join(
        settings.DISCOURSE_HOST, 'c', str(category['category']['id'])
    )

    logger.info(
        '- discourse.create_category(%s) topic_url=%s', category_name, topic_url
    )
    return category['category']['id'], topic_url


def process_category(category_details):
    """Check category is present in Table.
    If not create it in Discourse and store the id returned.
    """

    # DISCOURSE_DEV_POST_SUFFIX is used to differentiate the same target name
    # from different dev systems in Discourse.
    # It is not intended to be used for production when there is a dedicated Discourse.
    category_name = (
        category_details['category_name'] + settings.DISCOURSE_DEV_POST_SUFFIX
    )

    try:
        category = DiscourseCategory.objects.get(category_name=category_name)
        category_id = category.discourse_category_id
        post_url = os.path.join(settings.DISCOURSE_HOST, 'c', str(category_id))
    except DiscourseCategory.DoesNotExist:
        # Create Discourse client for fragalysis user
        client = _get_client()
        category_id, post_url = create_category(
            client,
            category_name,
            category_details['parent_name'],
            category_details['category_colour'],
            category_details['category_text_colour'],
        )
        DiscourseCategory.objects.create(
            category_name=category_name, discourse_category_id=category_id
        )
    return category_id, post_url


def create_post(user, post_details, category_id=None, topic_id=None):
    """Call discourse API posts.json to create a new Topic or Post within a topic"""
    assert user
    assert 'title' in post_details

    title = post_details['title']
    logger.info('+ discourse.create_post(%s, %s)', user.username, title)

    if not user.is_authenticated:
        return True, 'Please logon to Post to Discourse', 0, ''

    # Check user is present in Discourse
    client = _get_client()
    error, error_message, user_id = get_user(client, user.username)
    if user_id == 0:
        return error, error_message, 0, ''

    content = post_details['content']
    tags = post_details['tags'] if 'tags' in post_details else []

    if len(content) < 20:
        return True, 'Content must be more than 20 characters long in Discourse', 0, ''

    # Try to create a discourse topic post (a topic).
    # This might fail, especially if the topic title already exists.
    client = _get_user_client(user)
    try:
        post = client.create_post(content, category_id, topic_id, title, tags)
    except pydiscourse.exceptions.DiscourseClientError as dex:
        return True, dex.message, 0, ''

    # A topic's url is {URL}/t/{topic_id}/{post_number}
    post_url = os.path.join(
        settings.DISCOURSE_HOST, 't', str(post['topic_id']), str(post['post_number'])
    )

    logger.info(
        '- discourse.create_post(%s, %s) post_url=%s', user.username, title, post_url
    )
    return error, error_message, post['topic_id'], post_url


def process_post(category_id, post_details, user):
    """Check topic is present in Discourse.
    If it exists then post, otherwise create new topic for category
    """
    assert category_id
    assert 'title' in post_details

    error = False
    error_message = ''

    # DISCOURSE_DEV_POST_SUFFIX is used to differentiate the same target name
    # from different dev systems in Discourse.
    # It is not intended to be used for production when there is a dedicated Discourse.
    post_details['title'] = post_details['title'] + settings.DISCOURSE_DEV_POST_SUFFIX

    if topic := DiscourseTopic.objects.filter(
        topic_title=post_details['title']
    ).first():
        topic_id = topic.discourse_topic_id
        if post_details['content'] == '':
            # No content - Return the URL for the topic
            post_url = os.path.join(settings.DISCOURSE_HOST, 't', str(topic_id))
        else:
            # Create post for topic
            error, error_message, _, post_url = create_post(
                user, post_details, topic_id=topic_id
            )
    else:
        # Create Topic for Category
        error, error_message, topic_id, post_url = create_post(
            user, post_details, category_id=category_id
        )
        if not error:
            DiscourseTopic.objects.create(
                topic_title=post_details['title'],
                author=user,
                discourse_topic_id=topic_id,
            )

    return error, error_message, topic_id, post_url


def topic_posts(client, topic_id):
    """Gets posts for a given topic_id"""
    logger.info('+ discourse.topic_posts(%s)', topic_id)

    posts = client.topic_posts(topic_id)

    logger.info('- discourse.topic_posts(%s) posts=%s', topic_id, posts)
    return posts


def create_discourse_post(user, category_details=None, post_details=None):
    """Call discourse APIs to create the category, topic or post details as requested

    Makes the translation from Fragalysis topic_name to discourse id.
    :param user
    :param category_details (optional) (name=<target name required>, color=<hex> (optional defaults),
                                 text_color=<hex> (optional defaults), parent_category_name=<parent category>
                                 (initially can default to “Fragalysis targets” (this is a setting)
    :param post_details (optional) {title=<e.g session-project.name required>, content=<e.g session-project.description,
                                  created_date=<date>, tags[list]}
    :return: error or url of post.
    """
    assert user

    logger.info('+ discourse.create_discourse_post(%s)', user.username)

    error = False
    error_message = ''

    post_url = ''
    if category_details is None and post_details is None:
        # Nothing to post
        error_message = 'Please supply either category or post details'
        return True, post_url, error_message

    category_id = 0

    # If category details exist then create/return the category using the fragalysis user
    if category_details:
        category_id, post_url = process_category(category_details)

    # If post details exist then create the post
    if post_details:
        error, error_message, _, post_url = process_post(
            category_id, post_details, user
        )

    logger.info(
        '- discourse.create_discourse_post(%s) error=%s error_message="%s" post_url=%s',
        user.username,
        error,
        error_message,
        post_url,
    )
    return error, post_url, error_message


def list_discourse_posts_for_topic(topic_title):
    """Call discourse APIs to retrieve posts for the given topic

    Makes the translation from Fragalysis topic_name to discourse id.

    :param topic_title
    :return: error or json of posts.
    """
    print('+ discourse.get_discourse_posts_for_topic(%s)', topic_title)
    posts = ''
    error = False

    # DISCOURSE_DEV_POST_SUFFIX is used to differentiate the same target name
    # from different dev systems in Discourse.
    # It is not intended to be used for production when there is a dedicated Discourse.
    if topic_title:
        topic_title = topic_title + settings.DISCOURSE_DEV_POST_SUFFIX

    if topic := DiscourseTopic.objects.filter(topic_title=topic_title).first():
        # Create Discourse client
        client = _get_client()
        posts = topic_posts(client, topic.discourse_topic_id)
    else:
        error = True

    print('- discourse.get_discourse_posts_for_topic(%s) error=%s', topic_title, error)
    return error, posts


def check_discourse_user(user):
    """Call discourse API to check discourse user exists.
    If the user does not exit user_id will be 0.
    """

    client = _get_client()
    error, error_message, user_id = get_user(client, user.username)
    return error, error_message, user_id

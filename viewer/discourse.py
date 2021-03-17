"""
discourse.py Functions for to send and retrieve posts to/fraom the discourse platform.
"""
import os
from pydiscourse import DiscourseClient
from django.conf import settings
from viewer.models import DiscourseCategory, DiscourseTopic
import logging
logger = logging.getLogger(__name__)

def get_user(client, username):
    """Call discourse API users to retreive user by username.
    """
    logger.info('+ discourse.get_user')
    error = False
    error_message = ''
    user = 0

    try:
        user = client.user(username)
    except Exception as e:
        # User is not known in Discourse or there is a failure accessing Discourse.
        logger.error('discourse.get_user', 'get_user', exc_info=e)
        error = True
        error_message = 'Error validating user in Discourse. If this is your first post please log on to ' \
                        'Discourse once to create a User. URL is: ' + settings.DISCOURSE_HOST
    else:
        user = user['id']

    logger.info('- discourse.get_user')
    return error, error_message, user


def create_category(client, category_name, parent_name, category_colour='0088CC', category_text_colour='FFFFFF'):
    """Call discourse API categories.json to create a new (sub)category.
    """
    logger.info('+ discourse.create_category')
    category = client.create_category(category_name, color=category_colour,
                                      text_color=category_text_colour, parent=parent_name)
    topic_url = os.path.join(settings.DISCOURSE_HOST, 'c',  str(category['category']['id']))
    logger.info('- discourse.create_category')

    return category['category']['id'], topic_url


def process_category (category_details):
    """Check category is present in Table - If not create it in Discourse and store the id returned..
    """

    # DISCOURSE_DEV_POST_SUFFIX is used to differentiate the same target name from different dev systems in Discourse
    # It is not intended to be used for production when there is a dedicated Discourse.
    category_name = category_details['category_name'] + settings.DISCOURSE_DEV_POST_SUFFIX

    # Create Discourse client for fragalysis user
    client = DiscourseClient(
        settings.DISCOURSE_HOST,
        api_username=settings.DISCOURSE_USER,
        api_key=settings.DISCOURSE_API_KEY)

    try:
        category = DiscourseCategory.objects.get(category_name=category_name)
        category_id = category.discourse_category_id
        post_url = os.path.join(settings.DISCOURSE_HOST, 'c', str(category_id))
    except DiscourseCategory.DoesNotExist:
        category_id, post_url = create_category(client, category_name,
                                                category_details['parent_name'],
                                                category_details['category_colour'],
                                                category_details['category_text_colour'])
        DiscourseCategory.objects.create(category_name=category_name,
                                         discourse_category_id=category_id)
    return category_id, post_url


def create_post(user, post_details, category_id=None, topic_id=None):
    """Call discourse API posts.json to create a new Topic or Post within a topic
    """

    logger.info('+ discourse.create_post')
    if not user.is_authenticated:
        return True, 'Please logon to Post to Discourse', 0, ''

    # Create Discourse client for user
    client = DiscourseClient(
        settings.DISCOURSE_HOST,
        api_username=user.username,
        api_key=settings.DISCOURSE_API_KEY)

    # Check user is present in Discourse
    error, error_message, user_id = get_user(client, user.username)
    if user_id == 0:
        return error, error_message, user_id

    title = post_details['title']
    content = post_details['content']
    tags = post_details['tags']

    if tags is None:
        tags = []

    if len(content) < 20:
        return True, 'Content must be more than 20 characters long in Discourse', 0, ''

    post = client.create_post(content, category_id, topic_id, title, tags)
    # posts url = / t / {topic_id} / {post_number}
    post_url = os.path.join(settings.DISCOURSE_HOST, 't',  str(post['topic_id']),  str(post['post_number']))
    logger.info('- discourse.create_post')

    return error, error_message, post['topic_id'], post_url


def process_post(category_id, post_details, user):
    """Check topic is present in Discourse. IF exists then post, otherwise create new topic for category
    """
    error = False
    error_message = ''

    # DISCOURSE_DEV_POST_SUFFIX is used to differentiate the same target name from different dev systems in Discourse
    # It is not intended to be used for production when there is a dedicated Discourse.
    post_details['title'] = post_details['title'] + settings.DISCOURSE_DEV_POST_SUFFIX

    try:
        topic = DiscourseTopic.objects.get(topic_title=post_details['title'])
        topic_id = topic.discourse_topic_id
        if post_details['content'] == '':
            # No content - Return the URL for the topic
            post_url = os.path.join(settings.DISCOURSE_HOST, 't', str(topic_id))
        else:
            # Create post for topic
            error, error_message, null_id, post_url = create_post(user, post_details, topic_id=topic_id)
    except DiscourseTopic.DoesNotExist:
        # Create Topic for Category
        error, error_message, topic_id, post_url = create_post(user, post_details, category_id=category_id)
        if not error:
            DiscourseTopic.objects.create(topic_title=post_details['title'],
                                          author=user,
                                          discourse_topic_id=topic_id)
    return error, error_message, topic_id, post_url


def topic_posts(client, topic_id):
    """Gets posts for a given topic_id
    """
    logger.info('+ discourse.topic_posts')

    posts = client.topic_posts(topic_id)
    logger.info(posts)
    logger.info('- discourse.topic_posts')
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

    logger.info('+ discourse.create_discourse_post')
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
        error, error_message, topic_id, post_url = process_post(category_id, post_details, user)

    logger.info('- discourse.create_discourse_post')
    return error, post_url, error_message


def list_discourse_posts_for_topic(topic_title):
    """Call discourse APIs to retreive posts for the given topic

    Makes the translation from Fragalysis topic_name to discourse id.

    :param topic_title
    :return: error or json of posts.
    """
    print('+ discourse.get_discourse_posts_for_topic')
    posts = ''
    error = False

    # Get topic_id for title
    try:
        topic = DiscourseTopic.objects.get(topic_title=topic_title)
    except DiscourseTopic.DoesNotExist:
        topic = None

    if topic:
        # Create Discourse client
        client = DiscourseClient(
            settings.DISCOURSE_HOST,
            api_username=settings.DISCOURSE_USER,
            api_key=settings.DISCOURSE_API_KEY)
        posts = topic_posts(client, topic.discourse_topic_id)
    else:
        error = True

    print('- discourse.get_discourse_posts_for_topic')
    return error, posts

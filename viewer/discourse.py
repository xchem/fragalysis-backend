"""
discourse.py Functions for to send and retrieve posts to/fraom the discourse platform.
"""
import os
from pydiscourse import DiscourseClient
from django.conf import settings


def get_user(username):
    """Call discourse API users to retreive user by username.

    :param username
    :return: status_code and discourse user_id.
    """
    print('+ discourse.get_user')

    client = DiscourseClient(
        settings.DISCOURSE_HOST,
        api_username=settings.DISCOURSE_HOST,
        api_key=settings.DISCOURSE_API_KEY)

    user = client.user(username)
    print(user)
    print('- discourse.get_user')
    return user['id']


def create_category(category_name, parent_name, category_colour='0088CC', category_text_colour='FFFFFF'):
    """Call discourse API categories.json to create a new (sub)category.
    """

    print('+ discourse.create_category')
    client = DiscourseClient(
        settings.DISCOURSE_HOST,
        api_username=settings.DISCOURSE_HOST,
        api_key=settings.DISCOURSE_API_KEY)

    category = client.create_category(category_name, category_colour,
                                      text_color=category_text_colour, parent=parent_name)
    print(category)
    topic_url = os.path.join(settings.DISCOURSE_HOST, category['category']['topic_url'])
    print('- discourse.create_category')

    return category['category']['id'], topic_url


def create_post(content, tags=None, title=None, category_id=None, topic_id=None):
    """Call discourse API posts.json to create a new Topic or Post within a topic
    """
    if tags is None:
        tags = []

    print('+ discourse.create_post')

    client = DiscourseClient(
        settings.DISCOURSE_HOST,
        api_username=settings.DISCOURSE_HOST,
        api_key=settings.DISCOURSE_API_KEY)

    post = client.create_post(content, category_id, topic_id, title, tags)
    print(post)
    # posts url = / t / {topic_id} / {post_number}
    post_url = os.path.join(settings.DISCOURSE_HOST, 't',  str(post['topic_id']),  str(post['post_number']))
    print('- discourse.create_post')

    return post['topic_id'], post['post_number'], post_url


def topic_posts(topic_id):
    """Gets posts for a given topic_id
    """
    print('+ discourse.topic_posts')

    client = DiscourseClient(
        settings.DISCOURSE_HOST,
        api_username=settings.DISCOURSE_HOST,
        api_key=settings.DISCOURSE_API_KEY)

    posts = client.topic_posts(topic_id)
    print(posts)
    print('- discourse.topic_posts')
    return posts


def post_discourse_posts(category_details=[], post_details=[]):
    """Call discourse APIs to create the category, topic or post details as requested

    Makes the translation from Fragalysis topic_name to discourse id.

    :param category_details (optional) (name=<target name required>, color=<hex> (optional defaults),
                                 text_color=<hex> (optional defaults), parent_category=<parent category>
                                 (initially can default to “Fragalysis targets”)
    :param post_details (optional) (title=<e.g session-project.name required>, raw=<e.g session-project.description,
                                  created_date=<date>(defaults to current date))
    :return: error or url of post.
    """
    error = False
    url = ''

    return error, url


def get_discourse_posts_for_topic(topic_title):
    """Call discourse APIs to retreive posts for the given topic

    Makes the translation from Fragalysis topic_name to discourse id.

    :param topic_title
    :return: error or json of posts.
    """
    error = False
    posts = ''

    return error, posts

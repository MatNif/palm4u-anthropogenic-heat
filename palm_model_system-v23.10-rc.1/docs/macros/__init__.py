from verspec.python import Version
import os

from . LES_Model import add_palm_macros


docs_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def include_markdown(path):
    try:
        with open(os.path.join(docs_dir, 'content', path)) as f:
            content = f.read()
    except:
        content = '!!! warning\n    include_markdown(\''+path+'\') failed!\n'
    return content


def define_env(env):
    #chatter = env.start_chatting("define_env")
    #chatter('os.getcwd() = {}'.format(os.getcwd()))
    #print(os.getcwd())

    if 'mike' not in env.conf['plugins']:
        #chatter(f'extending site_url with "/local"')
        env.conf['site_url'] = env.conf['site_url'] + '/local'
    add_palm_macros(env)
    env.macro(include_markdown)

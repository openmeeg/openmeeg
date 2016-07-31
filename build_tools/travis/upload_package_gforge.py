import os
import sys
import glob
import paramiko

server = 'scm.gforge.inria.fr'
username = 'travisci'
password = os.environ['deploy_password']
remotepath = '/home/groups/openmeeg/htdocs/download/'
key_filename = os.path.expanduser('~/.ssh/openmeeg_deploy_key')

if __name__ == '__main__':
    ssh = paramiko.SSHClient()
    ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(server, username=username, password=password, key_filename=key_filename)
    sftp = ssh.open_sftp()
    for fname in glob.glob(sys.argv[1]):
        sftp.put(fname, remotepath + fname)
    sftp.close()
    ssh.close()

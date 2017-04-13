import os
import sys
import glob
import paramiko
from OpenSSL.crypto import load_privatekey, FILETYPE_PEM, dump_privatekey

# PLEASE BE SURE THESE VARIABLES ARE CORRECTLY DEFINED
# variables
server = os.environ['deploy_host']
username = os.environ['deploy_user']
password = os.environ['deploy_password']
remotepath = os.environ['deploy_folder']
key_filename = os.path.join(os.environ['src_dir'],'build_tools','openmeeg_deploy_key.pem')

if __name__ == '__main__':
    if len(sys.argv)==1:
        print('nothing to upload')
        sys.exit(0)
    # load the encrypted file
    with open(key_filename,'rb') as f:
        priv = f.read()

    pkey = load_privatekey(FILETYPE_PEM, priv, passphrase=password)
    priv = dump_privatekey(FILETYPE_PEM, pkey)

    # hack on the header
    priv = priv.replace('BEGIN ','BEGIN RSA ')

    key_filename += '.dec'

    with open(key_filename,'w') as f:
        f.write(priv)

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(server, username=username, key_filename=key_filename)
    sftp = ssh.open_sftp()

    for i in sys.argv[1:]:
        for fname in glob.glob(i):
            fname2 = fname.replace('build\\','')
            print('uploading file: ' + fname + ' ...')
            sftp.put(fname, os.path.join(remotepath, fname2))

    sftp.close()
    ssh.close()
    os.remove(key_filename)

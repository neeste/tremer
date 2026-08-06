import ftplib
import os

password = None
if os.path.exists(".deploy_secret"):
    with open(".deploy_secret", "r") as f:
        password = f.read().strip()

ftp = ftplib.FTP("bonkachen.com")
ftp.login("tremer_deploy@bonkachen.com", password)
print("Current dir:", ftp.pwd())
print("Files:")
ftp.dir()

import os
import ftplib
import getpass

def get_local_files(directories):
    files_to_upload = []
    for d in directories:
        for root, dirs, files in os.walk(d):
            for f in files:
                if f == '.DS_Store': continue
                local_path = os.path.join(root, f)
                remote_path = local_path.replace(os.sep, '/')
                files_to_upload.append((local_path, remote_path))
    return files_to_upload

def deploy(host, user, password, directories):
    print(f"Connecting to {host} via FTP...")
    try:
        # Try secure FTP over TLS first
        ftp = ftplib.FTP_TLS(host)
        ftp.login(user, password)
        ftp.prot_p()
    except Exception as e:
        print(f"FTPS failed ({e}), falling back to standard FTP...")
        ftp = ftplib.FTP(host)
        ftp.login(user, password)
        
    print("Connected successfully!")
    
    files = get_local_files(directories)
    total = len(files)
    uploaded = 0
    skipped = 0
    
    for i, (local, remote) in enumerate(files):
        # Create remote directories if they don't exist
        remote_dir = os.path.dirname(remote)
        if remote_dir:
            dirs = remote_dir.split('/')
            current = ""
            for d in dirs:
                current = f"{current}/{d}" if current else d
                try:
                    ftp.mkd(current)
                except ftplib.error_perm:
                    pass # Directory already exists
        
        # Check size to skip uploading unchanged files
        try:
            remote_size = ftp.size(remote)
            local_size = os.path.getsize(local)
            if remote_size == local_size:
                skipped += 1
                continue
        except Exception:
            pass # File doesn't exist remotely or SIZE not supported
            
        print(f"[{i+1}/{total}] Uploading {remote}...")
        with open(local, 'rb') as f:
            ftp.storbinary(f'STOR {remote}', f)
        uploaded += 1
        
    ftp.quit()
    print(f"\nDeployment complete! Uploaded {uploaded} files, Skipped {skipped} unchanged files.")

if __name__ == "__main__":
    host = "bonkachen.com"
    user = "tremer_deploy@bonkachen.com"
    print(f"Deployment Target: {host}")
    print(f"FTP Username:      {user}")
    
    password = getpass.getpass("FTP Password: ")
    directories = ['admin', 'Neely_project', 'wasm']
    
    try:
        deploy(host, user, password, directories)
    except ftplib.error_perm as e:
        print(f"\nAuthentication Failed: {e}")
        print("Please check that your FTP password is correct.")
    except Exception as e:
        print(f"\nDeployment Error: {e}")

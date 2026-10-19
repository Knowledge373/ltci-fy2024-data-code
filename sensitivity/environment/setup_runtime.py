#!/usr/bin/env python3
"""Fetch SHA256-pinned official packages and build workspace-local R 4.6.1.
Linux Debian 13 amd64 only. Does not install OS packages or write to /usr or /etc.
"""
import argparse,datetime,hashlib,json,os,platform,shutil,subprocess,urllib.error,urllib.request
from pathlib import Path
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True,help='New directory for downloaded packages and expanded runtime');p.add_argument('--cache-root',type=Path,help='Optional existing download tree; every cached artifact is SHA256-verified before copying');p.add_argument('--prepare-only',action='store_true',help='Explicit partial mode: download, verify, and extract R/compiler; do not compile R packages');a=p.parse_args();out=a.out.resolve();assert not out.exists() or not any(out.iterdir()),'Use a new or empty runtime directory';assert platform.system()=='Linux' and platform.machine() in ('x86_64','amd64');osrelease=Path('/etc/os-release').read_text();assert 'ID=debian' in osrelease and 'VERSION_ID="13"' in osrelease,'Pinned binary runtime requires Debian 13 amd64; no fallback version is selected';out.mkdir(parents=True,exist_ok=True)
 for cmd in ['dpkg-deb','gcc','g++','make','sha256sum']:
  assert shutil.which(cmd),f'Missing prerequisite {cmd}; install the official Debian 13 build prerequisite before setup'
 lock=json.loads((HERE/'DOWNLOAD_LOCK.json').read_text());receipt={'started_utc':now(),'status':'RUNNING','download_lock_sha256':sha(HERE/'DOWNLOAD_LOCK.json'),'downloads':[]}
 def save(): (out/'SETUP_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
 save()
 try:
  for f in ['DOWNLOAD_LOCK.json','toolchain_downloads.json','prepare_runtime.py','run-r-local.sh','rebuild_from_preserved_sources.R','rebuild_runtime.sh','Makevars.low-memory']:
   shutil.copy2(HERE/f,out/f)
  shutil.copytree(HERE/'verification-source',out/'verification-source')
  for d in lock['downloads']:
   target=out/d['relative_path'];target.parent.mkdir(parents=True,exist_ok=True);cache=a.cache_root/d['relative_path'] if a.cache_root else None;source=d['official_source_url']
   if cache and cache.is_file():
    assert sha(cache)==d['sha256'],f'Cache hash mismatch: {cache}';shutil.copy2(cache,target);used='verified local cache'
   else:
    urls=[source]
    if '/src/contrib/' in source and '/Archive/' not in source:
     filename=source.rsplit('/',1)[1];pkg=filename.split('_',1)[0];urls.append(source.replace('/src/contrib/','/src/contrib/Archive/'+pkg+'/'))
    last=None
    for url in urls:
     try:
      with urllib.request.urlopen(url,timeout=120) as response,target.open('wb') as dest:shutil.copyfileobj(response,dest)
      used=url;last=None;break
     except urllib.error.HTTPError as e:
      last=e
      if e.code not in [404,410]:raise
    if last:raise last
   assert target.stat().st_size==d['bytes'] and sha(target)==d['sha256'],f'Download hash/size mismatch: {target}'
   receipt['downloads'].append({'path':d['relative_path'],'source':used,'sha256':d['sha256'],'verified':True});save()
  for f in ['run-r-local.sh','rebuild_runtime.sh']:(out/f).chmod(0o755)
  cmd=['python3',str(out/'prepare_runtime.py')] if a.prepare_only else [str(out/'rebuild_runtime.sh')]
  receipt['command']=cmd;save()
  with (out/'SETUP_CONSOLE.log').open('w') as log:subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
  receipt.update({'status':'PASS_PARTIAL_RUNTIME_EXTRACTED' if a.prepare_only else 'PASS_PINNED_RUNTIME_BUILT','completed_utc':now(),'launcher':str(out/'run-r-local.sh')});save();print(json.dumps(receipt,indent=2))
 except BaseException as e:
  receipt.update({'status':'FAILED','error':repr(e),'stopped_utc':now()});save();raise
if __name__=='__main__':main()

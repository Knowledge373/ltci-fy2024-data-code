#!/usr/bin/env python3
"""Prepare a relocatable workspace-only pinned R runtime from preserved official packages.
No installation to /usr or /etc. Supported runtime platform: Debian13 amd64.
"""
import hashlib,json,os,platform,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
assert platform.system()=='Linux' and platform.machine() in ('x86_64','amd64'), 'This preserved binary runtime requires Linux amd64; do not substitute versions.'
lock=json.loads((ROOT/'DOWNLOAD_LOCK.json').read_text())
for f in lock['downloads']:
 p=ROOT/f['relative_path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==f['sha256'],str(p)
subprocess.run(['dpkg-deb','-x',str(ROOT/'r-base-core.deb'),str(ROOT/'r-runtime')],check=True)
for d in json.loads((ROOT/'toolchain_downloads.json').read_text()):
 subprocess.run(['dpkg-deb','-x',str(ROOT/(d['Package']+'.deb')),str(ROOT/'toolchain')],check=True)
r=ROOT/'r-runtime';home=r/'usr/lib/R'
for name in ['Makeconf','Renviron','Renviron.site','Rprofile.site','ldpaths','repositories']:
 p=home/'etc'/name;source=r/'etc/R'/name if name!='Renviron' else home/'etc/Renviron.ucf'
 if p.is_symlink():p.unlink()
 p.write_text(source.read_text().replace('/usr/lib/R',str(home)).replace('/usr/share/R',str(r/'usr/share/R')))
p=home/'bin/R';s=p.read_text().replace('R_HOME_DIR="/usr/lib/R"',f'R_HOME_DIR="{home}"',1).replace('R_SHARE_DIR="/usr/share/R/share"',f'R_SHARE_DIR="{r}/usr/share/R/share"').replace('R_INCLUDE_DIR="/usr/share/R/include"',f'R_INCLUDE_DIR="{r}/usr/share/R/include"').replace('R_DOC_DIR="/usr/share/R/doc"',f'R_DOC_DIR="{r}/usr/share/R/doc"');p.write_text(s)
(ROOT/'r-library').mkdir(exist_ok=True);(ROOT/'toolchain/lib').mkdir(exist_ok=True)
for lib,target in [('blas','blas.so.3'),('lapack','lapack.so.3'),('gfortran','gfortran.so.5'),('quadmath','quadmath.so.0')]:
 p=ROOT/'toolchain/lib'/f'lib{lib}.so';p.unlink(missing_ok=True);p.symlink_to('/usr/lib/x86_64-linux-gnu/lib'+target)
p=ROOT/'toolchain/usr/bin/x86_64-linux-gnu-gfortran';p.unlink(missing_ok=True);p.symlink_to('x86_64-linux-gnu-gfortran-14')
print('Prepared exact official R4.6.1 and workspace-only compiler dependencies. Package library installation is separate.')

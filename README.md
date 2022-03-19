Ensure that the system is Linux based with xsmgenspec and xsmcomputeflex installed. Follow the following steps for installation of xsmdas:
## Installation of xsmdas software
#### Prerequisites for xsmdas
1. OS: Linux/Unix (CentOS 7.0+, Ubuntu 14.04+, RHEL 6.5+, Fedora 20.0+, SLED 11.0+,
OS X 10.13+)
2. Compiler: gcc 4.4+

#### Steps for xsmdas installation
1. Download the installation package of XSMDAS and CALDB from https://pradan.issdc.gov.
in/pradan/
2. Unzip the installation package ch2_xsmdas_yyyymmdd_vn.mm.zip to desired directory:
``` 
unzip ch2_xsmdas_yyyymmdd_vn.mm.zip
```
3. Set environmental variables by adding following lines to ~/.bashrc
```
export xsmdas= <path to xsmdas directory>/xsmdas
export PATH="$xsmdas/bin:$xsmdas/scripts:$PATH"
export LD_LIBRARY_PATH="$xsmdas/lib/":$LD_LIBRARY_PATH
export PFILES="$PFILES:$xsmdas/pfiles"
```
  here 'path to xsmdas directory' is to be replaced with the absolute path under which
  xsmdas directory is present
  source ~/.bashrc file as:
```
.~/.bashrc
```
If the user has installed any other package that uses PIL with PFILES environment defined
in /.bashrc, make sure that the other declaration are appending the pfile paths to the
environment variable.
If the user has HEASOFT installation, it is recommended that the command to source
the initialization script of HEASOFT, commonly aliased as heainit, if included in .bashrc
shall be after the above commands. This will create a local copy of Parameter files under
$HOME/pfiles.
4. Installation of libraries:
```
cd $xsmdas
./InstallLibs
```
This will compile cfitsio and pil libraries from source.
5. Installation of CALDB:
Download the caldb zip file provided along with XSMDAS. Unzip the package
ch2_xsm_caldb_yyyymmdd.zip to $xsmdas directory as:
```
unzip ch2_xsm_caldb_yyyymmdd.zip -d $xsmdas
```
The caldb files will be extracted to $xsmdas/caldb directory

6. Compilation:
Once the libraries are installed compile XSMDAS with
```
cd $xsmdas
make
```
7. Modify one of the scripts file:
Change line 1 of xsmcomputeflux.py to:
``` 
#!/ysr/bin/env python3.9
```
(Note: Most of the above part of Installation steps was taken from the official XSM Data Analysis Guide provided by Physics Research Laboratory Ahmedabad and can be accessed from: https://pradan.issdc.gov.in/pradan/protected/downloadFile/xsm/ch2_xsm_data_analysis_guide.pdf)
## Installation
### For Backend
Install Python dependencies 
``` bash
$ pip3 install -r requirements.txt
```

Start the development server by running 
``` bash
# dev server with hot reload at http://localhost:8080
$ python app.py
```
### For Frontend

``` bash
$ npm install
```

or

``` bash
$ yarn install
```

### Basic usage

``` bash
# dev server with hot reload at http://localhost:3000
$ npm start 

# if you use Node 17+ use this command instead of `npm start`
$ npm run start:n17 
```

or 

``` bash
# dev server with hot reload at http://localhost:3000
$ yarn start

# if you use Node 17+ use this command instead of `yarn start`
$ yarn start:n17 
```


Summary: GNU Scientific Library (GSL)
name: gsl
Packager: rosalia@lanl.gov
%define version 0.6
%define release 0
version: %{version}
release: %{release}
#Prereq: 
#requires: 
vendor: the GSL team
Distribution: research software
copyright: Copyright (C) 1997, 1998, 1999, 2000 the GSL team
source: gsl-%{version}.tar.gz
group: Libraries/Research
%define mybuildroot /var/tmp/%{name}-build
%define installroot /install-tmp
BuildRoot: %{mybuildroot}

%description
  GSL is a library for numerical analysis that aims to be complete and
to follow modern coding conventions, as well as lending itself to being
used in very high level languages (VHLLs).

%prep
%setup -c
echo "dude, mybuildroot is " %{mybuildroot}
echo "dude, installroot is " %{installroot}
echo "dude, RPM_BUILD_ROOT is " $RPM_BUILD_ROOT

%build
cd %{name}-%{version}; ./configure --prefix=/usr; make

%install
cd %{name}-%{version}; make install prefix=%{mybuildroot}/usr
#cd %{name}-%{version}; make install prefix=%{mybuildroot}/%{installroot}

%post

%postun

%files
%doc %{name}-%{version}/{NEWS,ChangeLog,KNOWN-PROBLEMS,MACHINES,README,AUTHORS,THANKS}
%doc /usr/info/gsl-ref*
#%doc /usr/info/gsl-ref.info
#%doc /usr/info/gsl-ref.info-1
#%doc /usr/info/gsl-ref.info-2
#%doc /usr/info/gsl-ref.info-3
#%doc /usr/info/gsl-ref.info-4
#%doc /usr/info/gsl-ref.info-5
#%doc /usr/info/gsl-ref.info-6
#%doc /usr/info/gsl-ref.info-7
#%doc /usr/info/gsl-ref.info-8
#%doc /usr/info/gsl-ref.info-9
#%doc /usr/info/gsl-ref.info-10
#%doc /usr/info/gsl-ref.info-11
#%doc /usr/info/gsl-ref.info-12
#%doc /usr/info/gsl-ref.info-13
#%doc /usr/info/gsl-ref.info-14
#%doc /usr/info/gsl-ref.info-15
#%doc /usr/info/gsl-ref.info-16
#%doc /usr/info/gsl-ref.info-17
/usr/bin/gsl-config
/usr/bin/gsl-histogram
/usr/bin/gsl-randist
/usr/lib/gsl
/usr/include/gsl

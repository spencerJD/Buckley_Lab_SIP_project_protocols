* **system info**
  * system summary
    * `uname [-m]`
    * `uname -a`
      * give OS info
  * Hardware info
    * `sudo dmidecode -q`
  * CPU info
    * `less /proc/cpuinfo`
  * memory usage
    * `free`
  * memory occupied by each user
    * `sudo smem -u`
  * table of memory by user
    * `sudo smem -u | sort -nrk4 | head -n10 | awk '{print $1,int($4/1024+0.5)}'`
    * Alternative: `ps hax -o rss,user | awk '{a[$2]+=$1;}END{for(i in a)print i" "int(a[i]/1024+0.5);}' | sort -rnk2`
  * memory info
    * `less /proc/meminfo`
  * memory for all jupyter notebooks  
    * `memory-by-pids $(pgrep -f jupyter | xargs)`    
    * `sudo ~/bin/memory-by-pids $(pgrep -f jupyter -l | grep nick | cut -f 1 -d ' ' | xargs)`  
  * hard drive partitioning
    * `lsblk`
  * Hard drive space
    * `df`
      * `du -sh /path/to/folder/`
    * `vmstat`
  * who is logged on
    * `who`
  * processes for a user
    * `ps -u [username]`
  * users on system
    * `users`
  * groups on system
    * `groups`
  * groups for a user
    * `groups the_user`
  * add **NEW** user to secondary group
    * `useradd -G {group-name} username`
  * add **EXISTING** user to secondary group
    * `usermod -a -G GROUPNAME USERNAME`
  * file with open permissions
    * `chmod 777 <filename>`
  * open permissions on all files in a directory
    * `sudo chmod -R 777 *`
  * giving all users read/write permissions for files  
    * `chmod a+rw`
  * permissions for whole group in folder
    * `chmod g+s`
  * inheriting permissions
    * `chmod g+s /folder/`
  * List PATH:
    * `echo $PATH`
  * List shell:
    * `echo $SHELL`
  * Change nice-ness of a running job (higher value, lower priority)
    * `renice -n 10 11111`

* **Find path of executable**
  * `which <command>`
  * `whereis <command>`

* **Root**
  * login as root
    * `sudo -s`
    * `su -`

* **ssh key pair (keygen)**
  * [git](http://inchoo.net/tools-frameworks/how-to-generate-ssh-keys-for-git-authorization/)
  * [other](https://help.github.com/articles/generating-ssh-keys)

* **Maintainence**
  * Rebooting
    * reboot now
      * `sudo reboot`		
    * schedule reboot for later
      * `shutdown -r time "message"`		
        * time can be '+\d' format (ie., +5 = 5 minutes) or in 'HH:MM' 
  * updating packages
    * `sudo apt-get update`
    * `sudo apt-get upgrade`

* **Add new user**
  * **Add user to server**
    * Add user under "username"
      * `sudo useradd <username>`
      * on doe, use: `sudo adduser <username>`
    * Set user password
      * `sudo passwd <username>`
  * **Add user to samba (for server mount with smb)**
    * `sudo smbpasswd -a <username>`
    
* **Installations**
  * **R**
    * For global use on Rstudio server
      * `sudo R`
      * `>install.packages('MY_PACKAGE', lib='/opt/R/3.3.2/lib/R/library')`
    * For personal repo
      * /usr/bin/R
      * >install.packages('MY_PACKAGE')
    * For global use in conda's version of R
      * if package can be found on anaconda cloud `sudo /opt/anaconda2/bin/conda install -c r MY_PACKAGE`
      * if not `sudo /opt/anaconda2/bin/R` then `>install.packages('MY_PACKAGE')`
    * For personal R library, install using command line R without sudo
      * `R` then `>install.packages('MY_PACKAGE')`
  * **Python**
    * `sudo /opt/anaconda2/bin/pip install MY_PACKAGE`
    * `sudo /opt/anaconda2/bin/conda install -c MY_CHANNEL MY_PACKAGE`
  * **Unix/Linux**
    * updating packages
      * `sudo apt-get update`
      * `sudo apt-get upgrade`
    * removing packages
      * `apt-get remove --purge package`
      * `apt-get clean`
    * refresh .profile
      * `. ~/.profile`
  * **3rd party software (executables)**
    * For single executables, place in: `/opt/bioinfo/bin/`
      * This directory should already be in the PATH variable
      * `echo $PATH` to check 
    * For directories of many executables (eg., a directory cloned from GitHub)
      * Move directory to `/opt/bioinfo/`
      * Add the directory to the PATH variable
        * For all users: `/etc/bash_pathaddons`
          * `sudo nano /etc/bash_pathaddons
          * Add line: `export PATH="$PATH:/opt/bioinfo/MY_NEW_DIRECTORY_CHANGE_THIS/"`

* **Notifications**
  * message users
    * `write`  
  * messages via email
    * `echo “body_of_email” | mailx -s “email_subject_line” email@domain.com`
    * OR `mailx -s “email_subject_line” email@somewhere.com < file.txt  #use text file; will be email body`
      * **WARNING:** for screen, use `sleep 100` as the last command to make sure that mailx finishes  

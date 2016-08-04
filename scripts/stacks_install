#-----------------------------
#Setting up and Installing Stacks:
#-------------------------

#Install LAMP following this website
https://www.digitalocean.com/community/tutorials/how-to-install-linux-apache-mysql-php-lamp-stack-on-ubuntu

#during MySQL installation select root password

#Install phpMyAdmin following this website
https://www.digitalocean.com/community/tutorials/how-to-install-and-secure-phpmyadmin-on-ubuntu-12-04

#Install additional required packages
sudo apt-get install php-mdb2 php-mdb2-driver-mysql 
sudo apt-get install libdbd-mysql-perl 
sudo apt-get install samtools 
sudo apt-get install libbam-dev 

#Install Stacks allowing stacks to read both BAM and FASTQ files
tar xfvz stacks_x.xx.tar.gz 
cd stacks_x.xx 
./configure --enable-bam --with-bam-include-path=/usr/include/samtools --with-bam-lib-path=/usr/lib
make
sudo make install

#run MySQL as root and create a user (or do this from within phpMyAdmin)
mysql -u root -p
mysql> GRANT ALL ON *.* TO 'stacks_user'@'localhost' IDENTIFIED BY 'stackspassword';
mysql> exit;

#edit the MySQL config file adding in the stacks_user and stackspassword
cd /usr/local/share/stacks/sql/
sudo cp mysql.cnf.dist mysql.cnf 
sudo nano mysql.cnf


#Enable Stacks interface

#create stacks.conf file
sudo nano /etc/apache2/stacks.conf

#populate file with this script
<Directory "/usr/local/share/stacks/php">
        Order deny,allow
        Deny from all
        Allow from all
        Require all granted
</Directory>

Alias /stacks "/usr/local/share/stacks/php"

#make the stacks.conf file visible to Apache2
sudo nano /etc/apache2/apache2.conf

#add the stacks.conf to the file in the Include section
Include /etc/apache2/stacks.conf

#Restart apache
sudo service apache2 restart

#edit PHP config file to allow access to the MySQL database
sudo cp /usr/local/share/stacks/php/constants.php.dist /usr/local/share/stacks/php/constants.php 

#edit the PHP config file adding in the stacks_user and stackspassword
sudo nano /usr/local/share/stacks/php/constants.php 

#now can access  the web interface via http://localhost/stacks

#--------------------------
#check this out
http://gbs-cloud-tutorial.readthedocs.org/en/latest/index.html


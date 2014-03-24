amplitude-marginaliser
======================

INSTALL INSTRUCTIONS
====================

To install the ampmarginaliser library first run:
>> autoreconf --force --install

Following this run:
>> ./configure --prefix=MY_INSTALL_PATH
where MY_INSTALL_PATH is the path where you want the library
to be installed. By default this will install to /usr/local.
Then run:
>> make
>> make install
to install in the given path.

To use this library when compiling add the following flags:
-I/MY_INSTALL_PATH/include -L/MY_INSTALL_PATH/lib -lampmarginaliser


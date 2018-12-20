
 Higher-Order Spectral Analysis (HOSA) Toolbox.
 Version 2.0.3 (R12 compliant)  27 Dec 2000


There are three steps to installing the HOSA toolbox.

 (a) Download the SETUP.EXE file to your machine.
 (b) Extract folders from SETUP.EXE
 (c) Set MATLAB path so that the HOSA toolbox is included in MATLAB'S
     search path, and the documentation is available from the Launch Pad.


(a) The SmartCert procedure will download the SETUP.EXE file 
    to your machine.  Remember the folder where you downloaded the file. 

(b) SETUP.EXE is a self-extracting file. You can extract the HOSA      
     toolbox folders in (at least) two different ways:

    1) Under Windows Explorer, or My computer:
       Navigate to the folder where you downloaded the file,
         then click on the file SETUP.exe

       To continue, See `Using WinZip' below 
    
       OR
 
    2) Click the Start button
       Click on `Run...'
       Click the 'Browse' button,
          Navigate to the folder where you downloaded the file
           then double-click on the file SETUP.EXE
       You are now back in the Run Tab;
           click the `OK' button 
    
    'Using WinZip' 
       
     The `WinZip Self-Extractor' window will pop up. 
     This will display a default folder, \hosa, where the files will be 
     	unzipped. 
     You can over-ride this default by clicking on the Browse button.
     A folder will be created if it does not exist. 
     Click on the `Unzip' button to proceed.
     Another `WinZip Self-Extractor' window should pop up with the message
        125 file(s) unzipped successfully.
     Click the `OK' button.
     Click the `Close' button.
     
     The HOSA toolbox files have now been installed in the folder
     you specified. There will be two sub-folders, HOSA and
     HOSADEMO.  The former contains all the production mfiles;
     the latter contains the mfiles for the demo, and the matfiles
     with data for the examples in the on-line manual. 
     In addition to the two folders, there will be three files

       README.TXT (this file)
       EULA.DOC   (end-user license agreeement)
       HOSA.PDF   (HOSA documentation)

   The file HOSA.PDF contains the complete HOSA Toolbox Manual.
   You can view and print the manual via a standard PDF browser,
   such as Adobe's Acrobat Reader.

c) Start Matlab Ver 6 (Release 12). 
   We will use Matlab's Path tool to add the HOSA folders to the
   Matlab search path. 
   
   In the Matlab window,  click on File, then on `Set Path ...' 
       Click the button ``Add with Subfolders ...''
       Navigate to the folder where you installed the HOSA toolbox
                    under item b) above
       Click on the HOSA toolbox folder and then click the OK button
       Click on `Save' 
       Click on `Close'
   
   To make sure that the toolbox has been successfully installed,
   type `hosaver' (without the quotes) in the Matlab command window,
   and hit the enter key. 
    The following text should be displayed: 

      Higher-Order Spectral Analysis Toolbox.
      Version 2.0.3 (R12 compliant)  27 Dec 2000

   `hosahelp'  will give one-line syntax for all the toolbox mfiles.
   `help hosa' will give a functional descrption of the toolbox mfiles.
   `hosademo'  will run the HOSA toolbox demo. 


   Close your Matlab session.
   Start Matlab again.
   You should now see  ``Higher-Order Spectral Analysis Toolbox'' listed in
       the Launch Pad.  Click on it, and then on `Documentation' to get
        on-line help.
   



       


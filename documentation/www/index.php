<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <meta http-equiv="Content-type" content="text/html; charset=utf-8" />
        <meta http-equiv="Content-Language" content="en-us" />

        <title>OpenMEEG Project - Athena - INRIA</title>

        <meta name="ROBOTS" content="ALL" />
        <meta name="Copyright" content="(c) 2007 Copyright content: Project Athena - INRIA Copyright design: Alexandre Gramfort" />

        <link href="style.css" rel="stylesheet" type="text/css" media="all" />

        <!-- import the DOM logic from external javascript files -->
        <!--<script type="text/javascript" src="js/util.js"></script>-->
        <script type="text/javascript" src="js/jquery.js"></script>
        <script type="text/javascript" src="js/interface.js"></script>
        <script type="text/javascript" src="js/floatting_window.js"></script>
        <script type="text/javascript" src="js/util.js"></script>


</head>

    <body>
 <LINK REL="SHORTCUT ICON" href="img/OMlogo_tiny_cropped2.gif">
   <!-- begin div.full -->
    <div id="full">

        <!-- begin div.main -->
        <div id="main">

            <div id="top_right">
                <a href="http://www-sop.inria.fr/athena/">Athena Project Home Page</a><br/>
                <script type="text/javascript" charset="utf-8">
                // <![CDATA[
                echo_email();
                // ]]>
                </script>

                <div id="menu">
                    <a href="index.php?page=home">Home </a><br/>
                    <a href="https://gforge.inria.fr/frs/?group_id=435">Downloads (from forge)</a><br/>
                    <a href="index.php?page=why">Why use OpenMEEG?</a><br/>
                    <a href="index.php?page=documentation">Documentation</a><br/>
                    <a href="index.php?page=guide">Hacker's guide</a><br/>
                    <a href="index.php?page=publications">Publications</a><br/>
                    <a href="index.php?page=support">Support</a><br/>
                </div> <!-- end div.menu -->
            </div> <!-- end div.top_right -->

            <div>
                <div style="float:left">
                    <a href="."><img src="img/OMlogo_small.gif" width="150" height="145" alt="OMlogo Small">
                    </a>
                </div>
                <div id="title">
                    OpenMEEG
                    <div id="subtitle">
       A C++ package for low-frequency bio-electromagnetism <br>
     solving forward problems in the field of EEG and MEG.
                    </div> <!-- end div.subtitle -->
                </div>
            </div> <!-- end div.title -->


            <?php
                $pageOK = array(    'dev' => 'dev.php',
                                    'doc' => 'doc.php',
                                    'download' => 'download.php',
                                    'home' => 'home.php',
                                    'documentation' => 'documentation.php',
                                    'support' => 'support.php',
                                    'features' => 'features.php',
                                    'morenews' => 'morenews.php',
                                    'phpinfo' => 'phpinfo.php',
                                    'guide' => 'guide.php',
                                    'publications' => 'publications.php',
                                    'why' => 'why.php'
                               );

                if ( (isset($_GET['page'])) && (isset($pageOK[$_GET['page']])) ) {
                    include($pageOK[$_GET['page']]);   // Nous appelons le contenu central de la page
                }
                else {
                    include('home.php');   // Page par défaut quant elle n'existe pas dans le tableau
                }

            ?>

        </div>
        <!-- end div.main -->

        <!-- begin div.footer -->
        <div id="footer">
            <a href="http://validator.w3.org/check/referer">
                <img src="img/xhtml_valid.png" alt="valid xhtml" class="" />
            </a>
            <a href="http://jigsaw.w3.org/css-validator/check/referer">
                <img src="img/css_valid.png" alt="css xhtml" class="" />
            </a>
        </div>
        <!-- end div.footer -->

    </div>
    <!-- end div.full -->

    </body>
</html>




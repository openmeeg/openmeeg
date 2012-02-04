<div class="selected_heading">
    Home
</div>

<div class="part_title">
    What is OpenMEEG?
</div>

<div class="part_content" style="text-align:justify">
    <ul>
        <li>
           The <b>OpenMEEG</b> software is developed within the <a href="http://www-sop.inria.fr/athena" title="Athena Project Team">Athena project-team</a> at INRIA Sophia-Antipolis.
        </li>
        <li>
            In was initiated in 2006 by the Odyssee Project Team (INRIA/ENPC/ENS Ulm).
        </li>
        <li>
         OpenMEEG solves <b>forward problems</b> related to Magneto- and Electro-encephalography (MEG and EEG).
        </li>
        <li>
        The forward problem uses the symmetric Boundary Element Method (<strong>symmetric BEM</strong>), providing <a href="index.php?page=why">excellent accuracy</a>. (see <a href="index.php?page=publications">publications list</a>)
        </li>
        <!-- <li>
        The inverse problem  focusses on distributed source models with several types of regularization : the Minimum Norm (<strong>MN</strong>), the L2 of the surface gradient (<strong>HEAT</strong>) and the L1 of the surface gradient of the sources (<strong>TV</strong>).
        </li> -->
    </ul>
</div>

<div class="part_title">
    Using OpenMEEG
</div>
<div class="part_content">

    <div style="float:right">
        <a href="http://www.python.org/"><img src="img/python-logo.gif" width="211" height="71" alt="Python Logo"></a>
        <a href="http://fieldtrip.fcdonders.nl/" title="Fieldtrip"><img src="img/logo_fieldtrip.png" width="120" alt="Logo Fieldtrip"></a>
    </div>

    <div>
        <ul>
            <li>For M/EEG forward modeling:
                <ul>
                    <li>
                        From Matlab
                        for <a href="fieldtrip/openmeeg_eeg_leadfield_example.html">EEG</a>
                        and <a href="fieldtrip/openmeeg_meg_leadfield_example.html">MEG</a>
                        using <a href="http://fieldtrip.fcdonders.nl/" title="FieldTrip">Fieldtrip</a>.
                    </li>
                </ul>
            </li>
            
            <li>For general lead fields computation
                    (EEG, MEG, EIT, Internal potential):
                <ul>
                    <li>
                        From <a href="examples/python.html">Python</a>
                        (wrapping done using Swig)
                    </li>
                    <li>
                        From a <a href="examples/bash_script.html">Bash
                        script</a> on Unix systems.
                    </li>
                    <li>
                        From a <a href="examples/windows_bat.html">BAT file
                        </a>on Windows.
                    </li>
                </ul>
            </li>
            
        </ul>

    </div>

</div>

<div class="part_title">
    License
</div>
<div class="part_content">

    <div style="font-size:90%">
    OpenMEEG is distributed under the French opensource license <a
    href="http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html"
    title="CeCILL-B FREE SOFTWARE LICENSE AGREEMENT">CeCILL-B</a>. It is
    intended to give users the freedom to modify and redistribute the
    software. It is therefore compatible with popular opensource licenses
    such as the GPL and BSD licenses. The CeCILL-B license imposes to
    anybody distributing a software incorporating OpenMEEG the obligation
    to give credits (by citing the appropriate publications), in order for
    all contributions to be properly identified and acknowledged.
    </div>

    <p style="font-weight:bold;font-size:90%">The references to be acknowledged are:</p>
    <div class="publication">
        A. Gramfort, T. Papadopoulo, E. Olivi, M. Clerc. <strong>OpenMEEG: opensource software for quasistatic bioelectromagnetics,</strong>  <a href="http://www.biomedical-engineering-online.com/content/9/1/45" title="OpenMEEG: opensource software for quasistatic bioelectromagnetics"> BioMedical Engineering OnLine 45:9, 2010</a>
    </div>
    <div class="publication">
        Kybic J, Clerc M, Abboud T, Faugeras O, Keriven R, Papadopoulo T.
        <strong>A common formalism for the integral formulations of the forward EEG problem</strong>. <a  href=" http://ieeexplore.ieee.org/xpls/abs_all.jsp?isnumber=30034&arnumber=1375158&count=10&index=1">IEEE Transactions on Medical Imaging, 24:12-28, 2005.</a> <a  href="ftp://ftp-sop.inria.fr/odyssee/Publications/2005/kybic-clerc-etal:05.pdf">[PDF]</a>
    </div>

</div>


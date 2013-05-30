function echo_email ()
{
    var who_email = "openmeeg-devel"; 
    var who_domain = "lists.gforge.inria.fr";
    document.write("<a href=\"mailto"); 
    document.write(":" + who_email + "@" + who_domain +"\">"); 
    document.write("@ Contact us<\/a>");
}

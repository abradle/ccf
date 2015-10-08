
var ClassID = "ActiveIcm.ActiveIcmCtl";
var ICB_MIME_TYPE = "application/x-molsoft-icb";
var PLUGIN_KEYWORD_NAME = "molsoft";
var WIDTH  = "400";
var HEIGHT = "400";
var ACTIVEICM_LATEST_VERSION_URL = 'static/MMPMaker/activeICM.xml';
var INSERT_INTO = "body";

var httpRequest;
var os, browser;


var currentVersion = -1;
var icmPluginEnabled = false;
var icmPluginNotAccessible = true;
var xml_version;
var latestVersion;
var datapack_inserted = false;
var message_inserted = false;
var status_message_inserted = false;

function get_latest_version(callback1){
    if (window.XMLHttpRequest) { // Mozilla, Safari, ...
        httpRequest = new XMLHttpRequest();
        httpRequest.onreadystatechange = function(callback){
            when_ready(callback1);
        };
    } else if (window.ActiveXObject) { // IE
        httpRequest = new ActiveXObject("Microsoft.XMLHTTP");
        if (httpRequest) {
            httpRequest.onreadystatechange = function(callback){
                when_ready(callback1);
            };
        }
    }
    //asynchronous
    httpRequest.open('GET', ACTIVEICM_LATEST_VERSION_URL, true);
    httpRequest.send(null);
}

function process_response(responseXML, callback){//substitute     dump
    var x = responseXML.getElementsByTagName(xml_version);
    for (i = 0; i < x.length; i++) {
        var attributes = x[i].attributes;
        if ((attributes[1].value == os) && (attributes[0].value == browser)) {
            latestVersion = x[i].childNodes[0].nodeValue;
            if(callback){
                callback();

            }
        }
    }
}

function when_ready(callback){
    if (httpRequest.readyState == 4) {
        if(httpRequest.status == 200){
            process_response(httpRequest.responseXML, callback);
        }
        else
            window.alert(httpRequest.status);
    }
}

function detectBrowser() {
    var browser = "other";
    var userAgent = navigator.userAgent.toLowerCase();
    
    if(userAgent.indexOf("msie") != -1){
        browser = "msie";
    }
    else
    if (userAgent.indexOf("firefox") != -1) {
        browser = "firefox";
    }
    else
    if (userAgent.indexOf("chrome") != -1) {
        browser = "chrome";
    }
    else
    if (userAgent.indexOf("safari") != -1) {
        browser = "safari";
    }
    return browser;
}

function detectOS(){
    var os = "other";
    var userAgent = navigator.userAgent.toLowerCase();
    if(userAgent.indexOf("windows") != -1){
        os = "windows";
    }
    if(userAgent.indexOf("mac") != -1){
        os = "mac";
    }
    if(userAgent.indexOf("linux") != -1){
        os = "linux";
    }
    return os;
}


function getCurrentVersion(){    
    ppversion = -1;
    if( browser == "msie" ){
        document.write('<script type=\"text/vbscript\"> \n');
        document.write('on error resume next \n');
        document.write('Dim icbObject \n');
        document.write('set icbObject = CreateObject("' + ClassID + '") \n');
        document.write('ppversion = icbObject.pluginVersion \n');
        document.write('<\/script> \n');

    }else{
        icbObjectVersion = document.createElement("div");
        icbObjectVersion.setAttribute("style", "visibility:hidden");
        icbObjectVersion.innerHTML = "<OBJECT ID=\"icbObjectVersion\" type=\"" + ICB_MIME_TYPE+ "\" WIDTH=\"" + WIDTH + "px\" HEIGHT=\"" + HEIGHT + "px\" ></OBJECT>";
        document.body.appendChild(icbObjectVersion);
        ppversion = document.getElementById('icbObjectVersion').pluginVersion;
        if(detectBrowser()!="safari" || (detectOS()== "linux" && detectBrowser()!="firefox"))
            document.body.removeChild(icbObjectVersion);
    }
    return ppversion;
}
            
function isIcmPluginEnabled(){
    if(browser == "msie"){
        result = false;
        document.write('<script type=\"text/vbscript\"> \n');
        document.write('on error resume next \n');
        document.write('result = IsObject(CreateObject("' + ClassID + '")) \n');
        document.write('</SCRIPT> \n');
        return result;
    }else{
        if (navigator.plugins && navigator.plugins.length > 0) {
            for (i=0; i < navigator.plugins.length; i++ ) {
                if (navigator.plugins[i].name.toLowerCase().indexOf(PLUGIN_KEYWORD_NAME)!= -1){
                    return true;
                }
            }
        }
        return false;
    }
}


function insertDatapackObject(datapackFileName){
    //document.body.style.backgroundColor = "black";
    var ActiveIcmCtl = document.createElement("div");
    if(currentVersion <= 111){
        ActiveIcmCtl.innerHTML = "<OBJECT ID=\"ActiveIcmCtl\" type=\"" + ICB_MIME_TYPE+ "\" data=\"./" + datapackFileName + "\" WIDTH=\"100%\" HEIGHT=\"100%\" ></OBJECT>";
    }else{
        ActiveIcmCtl.innerHTML = "<OBJECT ID=\"ActiveIcmCtl\" type=\"" + ICB_MIME_TYPE+ "\" WIDTH=\"" + WIDTH + "px\" HEIGHT=\"" + HEIGHT + "px\" >" + "<param name=\"projectFile\" value=\"" +datapackFileName +"\" /></OBJECT>";
    }
    if(!datapack_inserted){//this is a problem with chrome and for some reason it runs this function twice
        if(INSERT_INTO!='body') {
		  document.getElementById(INSERT_INTO).appendChild(ActiveIcmCtl);
		}
		else {
		  document.body.appendChild(ActiveIcmCtl);
		}
        datapack_inserted = true;
    }
}

function insertDatapackObjectAfterCleanUp(datapackFileName){
    document.body.innerHTML = "";
    insertDatapackObject(datapackFileName);
}


function insertStatusBanner(){
    var statusDIV = document.createElement("div");
    statusDIV.setAttribute("name", "statusDIV");
    var innerHTML = "";
    if ( (browser != "other") && (os !="other")){
        innerHTML += " Operating system and Browser are supported.<br>";
    }else{
        innerHTML += " Operating system and Browser are NOT supported.<br>";
    }
    if (icmPluginEnabled){
        innerHTML +=  "activeICM plugin detected<br>";
    }else{
        innerHTML +=  "activeICM plugin NOT detected<br>";
    }
    if(icmPluginNotAccessible){
        innerHTML +=  "activeICM plugin not valid for the platform" + currentVersion + "<br>";
    }else{
        innerHTML +=  "activeICM plugin version " + currentVersion + "<br>";
    }
    innerHTML +=  "LATEST PLUGIN VERSION AVAILABLE " + latestVersion + "<br>";
    statusDIV.innerHTML = innerHTML;
    if(!status_message_inserted){
        document.body.appendChild(statusDIV);
        status_message_inserted = true;
    }
}

function insertIframe(src,id){
    /*http://stackoverflow.com/questions/1350767/accessing-javascript-variable-from-html-page-in-same-domain
     *http://www.w3schools.com/TAGS/tag_iframe.asp
     **/
    var iframe = document.createElement("iframe");
    iframe.setAttribute("src", src);
    iframe.setAttribute("id", id);
    iframe.setAttribute("frameborder",0);
    iframe.setAttribute("width","100%");
    iframe.setAttribute("height","100%");
    document.body.appendChild(iframe);
}

function append_message(id, contents, attributes){
    var message = document.createElement("div");
    message.setAttribute("id", id);
    if (attributes){
        for(i = 0; i< attributes.length;i++)
        message.setAttribute(attributes[i][0], attributes[i][1]);
    }
    message.innerHTML = contents;
    if (!message_inserted){
	    if(INSERT_INTO!='body') {
		  document.getElementById(INSERT_INTO).appendChild(message);
		}
		else {
          document.body.appendChild(message);
		}
        message_inserted = true;
    }
}


function compose_unsupported_platform_message(){
    return "<p align=\"center\"><b>Operating System or Browser not supported</b></p><table align=\"center\" bgcolor=\"#CC0000\" cellpadding =\"10\"><td><p>Contact us for further information: <a href=\"mailto:&#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107\"> &#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107</a></p></td></table>";

}

function compose_plugin_not_accessible_safari_mac_64_message(){
    return "<p align=\"center\"><b>Safari 4 in 64-bit mode detected</b></p><table bgcolor=\"#CC0000\" cellpadding =\"10\"><td><p>The activeICM plugin is already installed in your computer, but it is not yet supported for Safari 4 running in 64-bit mode</p><p>You can still use Safari 4 by running it in 32-bit mode:</p><p>Note: you'll need to quit Safari - print this page before continuing <input type=\"button\"; onClick=\"window.print()\"; value=\"Print\" /></p> <ol><li>Quit Safari (Cmd+Q)</li><li>On Finder, click on 'Go' and then 'Applications'</li><li>Find Safari and Ctrl-click on it</li><li>Pick 'Get Info'</li><li>Check the 'Open in 32-bit mode' box and close the dialog</li><li>Restart Safari 4</li></ol></p></td><td><a href=\"http://www.molsoft.com/getbrowser.cgi?product=activeicm&act=list\"><img style=\"border:none\" src=images/safari4-32bit.png width=\"200\"></a></td></table><p>After installation please restart  your browser.</p><p>Support: <a href=\"mailto:&#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107\"> &#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107 </a></p>";
}

function compose_new_plugin_IE(datapackFileName){
    var message = "<p align=\"center\"><b>Newer version of plugin available</b></p>";
    message += "<table bgcolor=\"#C7A317\" cellpadding =\"10\"><td>";
    message += "<p>You have an older version of the activeICM plugin installed in your computer.";
    message += "We recommend installing the latest version to ensure full functionality of this enhanced article.</p>";
    message += "<ul><li>Your version: "+currentVersion+"</li><li>Latest version: "+latestVersion+"</li></ul><p>";
    message += "Click on the button to the right to download the latest free activeICM plugin and to read the installation instructions.</p>";
    message += "  <p><b>Alternatively:</b> <a href=\"javascript:insertDatapackObjectAfterCleanUp('" + datapackFileName + "');\">Launch the enhanced article</a>";
    message += " using the present version of the plugin.</p></td><td><a  href=\"http://www.molsoft.com/getbrowser.cgi?product=activeicm&act=list\"><img style=\"border:none\" src=images/download_activeICM.png width=\"200\"></a></td></table><p>";
    message += "After installation please restart Internet Explorer. A secutiry bar will be displayed - please allow the activeICM control to run.</p><p>Support: <a href=\"mailto:&#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107\"> &#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107 </a></p>";
    return message;
}


function compose_new_plugin(datapackFileName){
    var message = "<p align=\"center\"><b>Newer version of plugin available</b></p><table bgcolor=\"#C7A317\" cellpadding =\"10\"><td>";
    message += "<p>You have an older version of the activeICM plugin installed in your computer.";
    message += "We recommend installing the latest version to ensure full functionality of this enhanced article.</p>";
    message += "<ul><li>Your version: "+currentVersion+"</li><li>Latest version: "+latestVersion+"</li></ul><p>";
    message += "Click on the button to the right to download the latest free activeICM plugin and to read the installation instructions.</p>";
    message += "  <p><b>Alternatively:</b> <a href=\"javascript:insertDatapackObjectAfterCleanUp('" + datapackFileName + "');\">Launch the enhanced article</a>";
    message += "using the present version of the plugin.</p></td>";
    message += "<td><a  href=\"http://www.molsoft.com/getbrowser.cgi?product=activeicm&act=list\"><img style=\"border:none\" src= \"images/download_activeICM.png\" width=\"200\"></a></td></table>";
    message += "<p>After installation please restart your browser.</p><p>Support: <a href=\"mailto:&#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107\"> &#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107 </a> </p>";
    return message;     
}

function compose_no_plugin_IE(){
    return "<p align=\"center\"><b>Visualisation plugin missing</b></p><table width=\"" + WIDTH + "\" bgcolor=\"#CC0000\" cellpadding =\"10\"><tr><td><a href=\"http://www.molsoft.com/getbrowser.cgi?product=activeicm&act=list\"><img style=\"border:none\" src=\"images/download_activeICM.png\" width=\"200px\"></a></td></tr><tr><td><p>A web plugin is required for the enhanced functionality of this article. Click on the button to the right to download the free activeICM plugin and to read the installation instructions.</p></td></tr></table><p>After installation please restart Internet Explorer. A security bar will be displayed - please allow the activeICM control to run.</p><p>Support: <a href=\"mailto:&#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107\"> &#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107 </a></p>";
}

function compose_no_plugin_safari_win(){
    return "<p align=\"center\"><b>Visualisation plugin missing</b></p><p align=\"center\"><b>***Click on CANCEL***</b><br>on the 'Safari can't find internet plug-in' window and follow the instructions below.</p><table width=\"" + WIDTH + "\" bgcolor=\"#CC0000\" cellpadding =\"10\"><tr><td><a href=\"http://www.molsoft.com/getbrowser.cgi?product=activeicm&act=list\"><img style=\"border:none\" src=images/download_activeICM.png width=\"200px\"></a></td></tr><tr><td><p>A web plugin is required for the enhanced functionality of this article. Click on the button to the right to download the free activeICM plugin and to read the installation instructions.</p></td></tr></table><p>After installation please restart  your browser.</p><p>Support: <a href=\"mailto:&#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107\"> &#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107 </a></p>";
}
function compose_no_plugin(){
    return "<p align=\"center\"><b>Visualisation plugin missing</b></p><table width=\"" + WIDTH + "\" bgcolor=\"#CC0000\" cellpadding =\"10px\"><tr><td><a href=\"http://www.molsoft.com/getbrowser.cgi?product=activeicm&act=list\"><img style=\"border:none\" src=\"images/download_activeICM.png\" width=\"200\"></a></td></tr><tr><td><p>A web plugin is required for the enhanced functionality of this article. Click on the button to the right to download the free activeICM plugin and to read the installation instructions.</p></td></tr></table><p>After installation please restart  your browser.</p><p>Support: <a href=\"mailto:&#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107\"> &#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107 </a></p>";
}
function compose_plugin_no_supported(){
    return "<p align=\"center\" style=\"color:black;background-color:white;\"><b>Visualisation plugin out-of-date</b></p><table width=\"" + WIDTH + "\" cellpadding =\"10px\" style=\"color:white;\"><tr><td><a href=\"http://www.molsoft.com/getbrowser.cgi?product=activeicm&act=list\"><img style=\"border:none\" src=\"images/download_activeICM.png\" width=\"200\"></a></td></tr><tr><td><p>Your version of activeICM plugin is no longer supported. Please update to the newest version available by clicking on the button to the right to download the free activeICM plugin and to read the installation instructions.</p></td></tr></table><p>After installation please restart  your browser.</p><p>Support: <a href=\"mailto:&#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107\"> &#105&#115&#101&#101&#64&#115&#103&#99&#46&#111&#120&#46&#97&#99&#46&#117&#107 </a></p>";
}

function generateICBDocument(datapackFileName){
    
    switch(true){
        case ((browser == "other") || (os =="other")):
            append_message("unsupported_platform", compose_unsupported_platform_message());
            break;
        case (icmPluginEnabled && icmPluginNotAccessible && browser =="safari" && os == "mac"):
            append_message("plugin_not_accessible",compose_plugin_not_accessible_safari_mac_64_message());
            break;
        case (icmPluginEnabled && icmPluginNotAccessible):
            append_message("no_plugin",compose_plugin_no_supported(),[["style","color:white; background-color:black;"]]);
            break;
        case (icmPluginEnabled && (latestVersion <= currentVersion)):
            insertDatapackObject(datapackFileName);
            break;
        case (icmPluginEnabled && (latestVersion > currentVersion) && (browser == "msie")):
            append_message("new_plugin_IE",compose_new_plugin_IE(datapackFileName));
            break;
        case (icmPluginEnabled && (latestVersion > currentVersion)):
            append_message("new_plugin",compose_new_plugin(datapackFileName));
            break;
        case (!icmPluginEnabled  && (browser == "msie")):
            append_message("no_plugin_IE",compose_no_plugin_IE());
            break;
        case (!icmPluginEnabled  && (browser == "safari" && os == "windows")):
            append_message("no_plugin_safari_win",compose_no_plugin_safari_win());
            break;
        default:
            append_message("no_plugin",compose_no_plugin());
    }
}



/*http://ajaxian.com/archives/serial-async-xhr*/
function composeICBdocument(datapackFileName,spec_width,spec_height,insert_into,xml_url) {
    WIDTH = spec_width;
    ACTIVEICM_LATEST_VERSION_URL = xml_url;
	HEIGHT = spec_height;
	INSERT_INTO = insert_into;
	getLocalParameters(function () {
        get_latest_version(function () {
            generateICBDocument(datapackFileName);
        });
    });
}
/*This function considers that the global variables browser, os and currentVersion have been assigned*/
function detect_pluginAccessibility(){
    icmPluginNotAccessible = icmPluginEnabled && (currentVersion == undefined || currentVersion == -1);
}

function getLocalParameters(callback2){
    browser = detectBrowser();
    os = detectOS();
    xml_version = (browser == "safari" || browser == "chrome") ?  "activeICM_latest_version" : "tns:activeICM_latest_version";
    icmPluginEnabled = isIcmPluginEnabled();
    currentVersion = -1;
    currentVersion = getCurrentVersion();
    detect_pluginAccessibility();
    if (callback2) {
        callback2();
    } 
}




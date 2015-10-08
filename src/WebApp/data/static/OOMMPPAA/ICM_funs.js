var mybase = location.protocol+"//"+location.hostname+":"+location.port;
if (location.port == ""){
  mybase = location.protocol+"//"+location.hostname;
};
// for displaying MMPs
var dispconst = 0;
function runCommand(s_command) {
var ctl =document.getElementById('ActiveIcmCtl');
ctl.RunCommands(s_command);
}

//using these jquery functions:
function getVar(varstr) {
      var ctl =this.document.getElementById('ActiveIcmCtl');
      return ctl.GetShellVar(varstr).replace(/\.,/ig,".0,");
}


function findasgr(s_me,mycanv){
if(typeof(mycanv)==='undefined'){
mycanv = 'mymolcanv';
}

runCommand('my_mol = as_graph &! a_'+target_name+'TEMP.*');
var my_mol = getVar('my_mol').replace(/ /g,'').replace(/[\n\r]/g, '');
var s_m = '{"my_mol":'+"[]"+'}'

runCommand('my_num = Nof(as_graph & (a_shape.* | a_act.* | a_inact.*))');
var my_num = getVar('my_num').replace(/ /g,'').replace(/[\n\r]/g, '');
var men = JSON.parse(my_num);
var my_num = men['my_num']
// Check to see if there's too many points (too many being 40)
if (my_num > 40){
  document.getElementById('MMPTooManyError').style.display='block'; 
}

else if (my_num == 0){
// here we have selected no spots
document.getElementById('MMPError').style.display='block'; 
}
else{
try{
document.getElementById('mymmpdesc').style.display='none';
document.getElementById('mymmpdesc2').style.display='block';
}
catch(err)
  {
  var meeeeee = 2;
  }

runCommand('my_res = String(Atom(as_graph))');
var my_var = getVar('my_res').replace(/ /g,'');

// put a loading gif into the MMP viewer
// load in the SDs and then load in the image
var url_string = mybase+'/Viewer/loader/?function=MOLS&choice='+my_var+'&map='+s_me+'&output=images';
url_string = url_string.replace(/[\r\n]/g, "");
//fillme(url_string,mycanv,new_s)
// Make the URL for the molecule loader
var url_string_new = mybase+'/Viewer/loader/?function=MOLS&choice='+my_var+'&map='+s_me+'&output=sds';
url_string_new = url_string_new.replace(/[\r\n]/g, "");
var new_s = "read mol '"+url_string_new+"'";
dispconst = 1;
fillme(url_string,mycanv,new_s)

}

}


// function to make a box

function makeBox(){
var my_var = getVar('as_graph').replace(/ /g,'');
// If it
var n = my_var.length;

if(n!=20)
{
runCommand('display box Box(as_graph)');
}
else{
alert("Nothing selected to make box. Please right-click drag to select atoms");
//
};
};

// function to find experiments potential MMPs with 

function findexps(s_me,mycanv){
if(typeof(mycanv)==='undefined'){
mycanv = 'mymolcanv';
}
runCommand('delete my_var');
runCommand('my_var = String(Box())');
var my_var = getVar('my_var').replace(/ /g,'');
var url_string = mybase+'/Viewer/loader/?function=EXP&choice='+my_var+'&map='+s_me+'&output=images';
url_string = url_string.replace(/[\r\n]/g, "");
fillme(url_string,mycanv)
dostuff();
}




function loadmymmp(mmpid,counter){
var myvar = counter.toString();
if ($('#choicecmp'+myvar).is(':checked')){
var target_id = '{{ target.id }}';
var url_string = mybase+'/Viewer/loader/?choice='+mmpid+'&function=GETMMP';

runCommand('delete a_.')
runCommand('read mol "'+url_string+'" name="s'+myvar+'"');
runCommand('display a_. xstick');
runCommand('display a_*_*.m//c* xstick purple');
runCommand('display a_.m//Na* cpk blue');
runCommand('display a_.m//Li* cpk red');
runCommand('display a_.m//k* cpk brown');
runCommand('display a_.m//Zn* cpk cyan');

}
else{
runCommand('delete a_s'+myvar+'.');
};
}



function dostuff(choice){
//runCommand('center static')

runCommand('GRAPHICS.l_redraw = no')

if (choice=="act"){
makeshape('a_act.na*','star','blue');
makeshape('a_act.li*','star','red');

makeshape('a_act.zn*','star','cyan');
makeshape('a_act.k*','star','brown');


};

if (choice=="inact"){
makeshape('a_inact.na*','cube','blue');
makeshape('a_inact.li*','cube','red');

makeshape('a_inact.zn*','cube','cyan');
makeshape('a_inact.k*','cube','brown');
}


if (choice=="non"){
makeshape('a_non.li*','doublepyramid','red');
makeshape('a_non.na*','doublepyramid','blue');
makeshape('a_non.zn*','doublepyramid','cyan');
makeshape('a_non.k*','doublepyramid','brown');
}

runCommand('GRAPHICS.l_redraw = yes')
runCommand('display new')
runCommand('display a_shape.na* purple xstick');
runCommand('color a_shape.na* Trim(Occupancy( a_shape.na*),0.0001,1.0001)//0.0001//1.0001'); 

}




function loadpdb(choice,pk){
var myvar = choice;
if ($('#choice'+myvar).is(':checked')){
var url_string = mybase+'/Viewer/loader/?choice='+pk+'&function=GETMAP';
runCommand('read pdb "'+url_string+'" name="'+choice+'"');

//runCommand('color background rgb={255,255,255}');
dostuff(choice)
}
else{
runCommand('delete a_'+choice+'.');
};

}



function loadpdblight(choice,pk){
var myvar = choice;

if ($('#choice'+myvar+'blah').is(':checked')){
//document.getElementById('nowwhat2').style.display='';
//document.getElementById('sliderdesc').style.display='block';
//document.getElementById('toggledesc').style.display='none';
var url_string = mybase+'/Viewer/loader/?choice='+pk+'&function=GETMAP';
runCommand('read pdb "'+url_string+'" name="'+choice+'"');
dostuff(choice)
//runCommand('center a_');
}
else{
if(choice=="act"){
	runCommand('deleteObjs "star" Atom(a_'+choice+'.)');
};
if(choice=="inact"){

	runCommand('deleteObjs "cube" Atom(a_'+choice+'.)');
};
if(choice=="non"){

	runCommand('deleteObjs "doublepyramid" Atom(a_'+choice+'.)');
};

runCommand('delete a_'+choice+'.');
};

}



function loadmyprot(pdb_code){
var url_string = mybase+'/Viewer/loader/?choice='+pdb_code+'&function=VIEWPROTEIN';
runCommand('read pdb "'+url_string+'" name="'+pdb_code+'"');
runCommand('display a_ xstick');
};

function loadprot(pdb_code,counter){
var myvar = counter.toString();
if ($('#choice'+myvar).is(':checked')){
var url_string = mybase+'/Viewer/loader/?choice='+pdb_code+'&function=VIEWPROTEIN';
runCommand('read pdb "'+url_string+'" name="'+pdb_code+'"');
runCommand('display a_ xstick');
}
else{
runCommand('delete a_'+pdb_code+'.');
};
}

function loadmol2(pdb_code,dist){
//alert(dist);
var my_pdb = pdb_code.substring(0, 4);
runCommand('delete a_.*')
var url_string = mybase+'/Viewer/loader/?choice='+pdb_code+'&function=VIEWMOL';
runCommand('read mol "'+url_string+'"');
runCommand('display a_'+my_pdb+'*.* xstick');

}


function loadmol(pdb_code,counter){

var my_pdb = pdb_code.substring(0, 4);
var myvar = counter.toString();
if ($('#choicemol'+myvar).is(':checked')){
var url_string = mybase+'/Viewer/loader/?choice='+pdb_code+'&function=VIEWMOL';
runCommand('read mol "'+url_string+'"');
runCommand('display a_'+my_pdb+'*.* xstick');

}
else{
runCommand('delete  a_'+my_pdb+'*.*')
}
}


function loadcanvmols(pdb_code){


var url_string = mybase+'/Viewer/loader/?choice='+pdb_code+'&function=VIEWMOL';
runCommand('read mol "'+url_string+'"');
runCommand('display a_'+pdb_code+'*.* xstick');
var url_string_new =mybase+'/Viewer/loader/?function=MAKESIM&map=1&choice='+pdb_code;
fillme(url_string_new,'mynewercanv');
var url_string_new =mybase+'/Viewer/loader/?function=MAKESIM&map=2&choice='+pdb_code;
fillme(url_string_new,'mynewercanv2');
var url_string_new =mybase+'/Viewer/loader/?function=MAKESIM&map=3&choice='+pdb_code;
fillme(url_string_new,'mynewercanv3');
var url_string_new =mybase+'/Viewer/loader/?function=MAKESIM&map=4&choice='+pdb_code;
fillme(url_string_new,'mynewercanv4');
document.getElementById('mymoldiv').style.display='none';
document.getElementById('molfigs').style.display='';

}

function loadimagemol(pdb_code,counter){

var my_pdb = pdb_code.substring(0, 4);
var myvar = counter.toString();

var me =  document.getElementById('choicemol'+myvar).getElementsByTagName("img")[0]
var descdiv = document.getElementById('choicemol'+myvar).getElementsByClassName("mydescription")[0]
var alt = me.alt;
if (me.alt=="off"){
me.alt = "on";
descdiv.style.background = "#555";
descdiv.style.color = "white";
var url_string = mybase+'/Viewer/loader/?choice='+pdb_code+'&function=VIEWMOL';
runCommand('read mol "'+url_string+'"');
runCommand('display a_'+my_pdb+'*.* xstick');


//loadcanvmols(pdb_code)
//var url_string_new =mybase+'/Viewer/loader/?function=MAKESIM&choice='+pdb_code;
//fillme(url_string_new);

}

else{
me.alt="off";
descdiv.style.background = "";
descdiv.style.color = "#666";
runCommand('delete  a_'+my_pdb+'*.*')
}


}



function loadfrags(pdb_code,counter){

var myvar = counter.toString();
if ($('#choicefrag'+myvar).is(':checked')){
var url_string = mybase+'/Viewer/loader/?choice='+pdb_code+'&function=VIEWFRAG';
runCommand('read mol "'+url_string+'"');
runCommand('display a_. xstick');
}
else{
runCommand('delete a_'+pdb_code+'_*.')
}
}

function loadallmols(){

var url_string = mybase+'/Viewer/loader/?choice='+target_id+'&function=VIEWALLMOLS';
if ($('#allmols').is(':checked')){
runCommand('read mol "'+url_string+'"');
runCommand('display a_.* xstick');
}
else{
runCommand('delete a_.*')
}
}

function dssurface(){
if ($('#dssurf').is(':checked')){
runCommand('display a_ skin');
}
else{
runCommand('delete a_ ribbon')
}
}

function whatisasgraph(){
runCommand('delete a_. &! as_graph')
}

function loadgroup(groupchoice,choicegr){
var url_string = mybase+'/Viewer/loader/?choice='+target_id+'&groupchoice='+groupchoice+'&function=VIEW_MOL_GROUP';
if ($('#'+choicegr).is(':checked')){
runCommand('read mol "'+url_string+'"');
runCommand('display a_.* xstick');
}
else{
runCommand('delete a_.*')
}
}


// From the original script


// Then show the atoms as shapes
function makeshape(atoms,shape,color){
// first of all check to see if one if these exists
	var my_acts = $('#act').slider('getValue')[0].parentElement.childNodes[1].childNodes[1].innerHTML.split(":");
	var my_act = parseFloat(my_acts[0]).toFixed(2);
	var max_act = parseFloat(my_acts[1]).toFixed(2);
	var my_points = $('#pharm').slider('getValue')[0].parentElement.childNodes[1].childNodes[1].innerHTML.split(":");
	var my_point = parseInt(my_points[0]);
	var max_point = parseInt(my_points[1]);
	var my_s = 'showshapes "'+shape+'" "'+color+'" Atom('+atoms+') '+my_act+' '+my_point+' '+max_act+' '+max_point;
	runCommand(my_s);
};

var opts = {
  lines: 13, // The number of lines to draw
  length: 20, // The length of each line
  width: 10, // The line thickness
  radius: 30, // The radius of the inner circle
  corners: 1, // Corner roundness (0..1)
  rotate: 0, // The rotation offset
  direction: 1, // 1: clockwise, -1: counterclockwise
  color: '#000', // #rgb or #rrggbb or array of colors
  speed: 1, // Rounds per second
  trail: 60, // Afterglow percentage
  shadow: false, // Whether to render a shadow
  hwaccel: false, // Whether to use hardware acceleration
  className: 'spinner', // The CSS class to assign to the spinner
  zIndex: 2e9, // The z-index (defaults to 2000000000)
  top: '50%', // Top position relative to parent
  left: '50%' // Left position relative to parent
};


function fillme(urlstring,divid,new_s){
// Make an  image
  var img = document.createElement("IMG");
  img.src = urlstring;
// Make the spinner
  var spinner = new Spinner(opts).spin(document.getElementById(divid));

  if (new_s != undefined){
	var inputs = document.getElementsByTagName("INPUT");
	for (var i = 0; i < inputs.length; i++) {
	        inputs[i].disabled = true;
			}
	document.getElementById('updatemapbut').setAttribute("class","disabled btn btn-sm btn-success");
	// If  the cancel button is pressed then cancel this function
	
	    img.onload =  function(){
	      $('#'+divid).html('');
	      // refresh the dots
	      dostuff();
	      img.style.width="100%";
	      runCommand("GRAPHICS.l_redraw = no")
	      runCommand("delete a_*MOL*.");
	      runCommand("GRAPHICS.l_redraw = yes")
	      runCommand("display new");
	      runCommand(new_s);
	      runCommand("display a_*MOL*. xstick");
	      runCommand("colorMols Mol(a_MOL*.*)");
	      $('#'+divid).html('');
	      document.getElementById(divid).appendChild(img);
	      //runCommand("center Mol(a_MOL*.*)");
	for (var i = 0; i < inputs.length; i++) {
	        inputs[i].disabled = false;
	}
	document.getElementById('updatemapbut').setAttribute("class","btn btn-sm btn-success");
	
	      }
   }
  else{
    img.onload =  function(){
    
      spinner.stop();
      img.style.width="100%";
      document.getElementById(divid).appendChild(img);
       };
     };
  };






function toghelp(){
var helps = document.getElementsByClassName('description');
if ($('#helptog').is(':checked')){
for (var i = helps.length - 1; i >= 0; i--)
{
  helps[i].style.cssText = 'display:none !important';
}
}
else{
for (var i = helps.length - 1; i >= 0; i--)
{
  helps[i].style.cssText = 'display:visible';
}
}

};



var browser = "";
var OSName="Unknown OS";
if (navigator.appVersion.indexOf("Win")!=-1) OSName="Windows";
if (navigator.appVersion.indexOf("Mac")!=-1) OSName="MacOS";
if (navigator.appVersion.indexOf("X11")!=-1) OSName="UNIX";
if (navigator.appVersion.indexOf("Linux")!=-1) OSName="Linux";


  function isIcmPluginEnabled(){
    
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

        if (navigator.plugins && navigator.plugins.length > 0) {
            for (i=0; i < navigator.plugins.length; i++ ) {
                if (navigator.plugins[i].name.toLowerCase().indexOf('molsoft')!= -1){
                    return true;
                }
            }
        return false;
    }
}
<?php
function sanitizeString($var){
   if(get_magic_quotes_gpc()) $var = stripslashes($var);
   $var=htmlentities($var);
   $var=strip_tags($var);
   return $var;
}

function valueToJsonValue($value, $type){

   if($value == "null" || $value == ""){
      return "null";
   }

   if($type == "str"){
      $typed_value = strval($value);
   }

   if($type == "int"){
      $typed_value = intval($value);
   }

   if($type == "float"){
      $typed_value = floatval($value);
   }


   return json_encode($typed_value);
}

function createTxtFilecontent($args)
{

/*
{
  "FILE_FORMAT": 1, 
  "FILE_HEADER_LINES": 4, 
  "X1_FACTOR": 0.001, 
  "Y1_FACTOR": 1e-06, 
  "YEAR_FABRICATION": null, 
  "SURFACE_SHAPE": "plane", 
  "FUNCTION": null, 
  "LENGTH": 1.36, 
  "WIDTH": 0.16, 
  "THICK": 0.05, 
  "LENGTH_OPTICAL": 1.2, 
  "SUBSTRATE": "silicon", 
  "COATING": null, 
  "FACILITY": "ESRF", 
  "INSTRUMENT": null, 
  "POLISHING": null, 
  "ENVIRONMENT": null, 
  "SCAN_DATE": null, 
  "PLOT_TITLE_X1": "x(mm)", 
  "PLOT_TITLE_Y1": "Slope (urad)", 
  "CALC_HEIGHT_RMS": null, 
  "CALC_HEIGHT_RMS_FACTOR": null, 
  "CALC_SLOPE_RMS": null, 
  "CALC_SLOPE_RMS_FACTOR": null, 
  "USER_EXAMPLE": "This is an example of user keyword"
}
*/


$fields = array(
array("FILE_FORMAT","int"),
array("FILE_HEADER_LINES","int"),
array("X1_FACTOR","float"),
array("Y1_FACTOR","float"),
array("YEAR_FABRICATION","str"),
array("SURFACE_SHAPE","str"),
array("FUNCTION","str"),
array("LENGTH","float"),
array("WIDTH","float"),
array("THICK","float"),
array("LENGTH_OPTICAL","float"),
array("SUBSTRATE","str"),
array("COATING","str"),
array("FACILITY","str"),
array("INSTRUMENT","str"),
array("POLISHING","str"),
array("ENVIRONMENT","str"),
array("SCAN_DATE","str"),
array("PLOT_TITLE_X1","str"),
array("PLOT_TITLE_Y1","str"),
array("CALC_HEIGHT_RMS","str"),
array("CALC_HEIGHT_RMS_FACTOR","str"),
array("CALC_SLOPE_RMS","str"),
array("CALC_SLOPE_RMS_FACTOR","str"),
array("USER_EXAMPLE","str")
);

$txt_content = "{\n";//array("{");

$number_fields = count($fields);

$cur_line = 0;
foreach ($fields as $line){
    $value = valueToJsonValue($args[$line[0]],$line[1]);
    $txt_line = "\"" . $line[0] . "\": " . $value;

    if($cur_line+1 < $number_fields){
       $txt_line = $txt_line . ",";
    } 
//    $txt_line = sanitizeString($txt_line);
    $txt_line = "  " . $txt_line;
//    array_push($txt_content, $txt_line);
    $txt_content .= $txt_line . "\n";
    $cur_line+=1;
}
    //array_push($txt_content, "}");
    $txt_content .= "}\n";

return $txt_content;
}

function createEmailContent($message, $datafile, $txt_content){

    // ordinary message $message
    // array with filenames to be sent as attachment $files

    // boundary 
    $semi_rand = md5(time()); 
    $mime_boundary = "==Multipart_Boundary_x{$semi_rand}x"; 

    // multipart boundary 
    $message = "This is a multi-part message in MIME format.\n\n" . "--{$mime_boundary}\n" . "Content-Type: text/plain; charset=\"iso-8859-1\"\n" . "Content-Transfer-Encoding: 7bit\n\n" . $message . "\n\n"; 
    $message .= "--{$mime_boundary}\n";

    // preparing attachments
    $datafile_name = "dabam.dat";
    $txtfile_name = "dabam.txt";

    // data file
    $file = fopen($datafile,"rb");
    $data = fread($file,filesize($datafile));
    fclose($file);
    $data = chunk_split(base64_encode($data));
    $message .= "Content-Type: text/plain; charset=UTF-8;\n" . " name=\"$datafile_name\"\n" . 
        "Content-Disposition: attachment;\n" . " filename=\"$datafile_name\"\n" . 
        "Content-Transfer-Encoding: base64\n\n" . $data . "\n\n";
    $message .= "--{$mime_boundary}\n";

    // txt file
    $data = $txt_content;
    $data = chunk_split(base64_encode($data));
    $message .= "Content-Type: text/plain; charset=UTF-8;\n" . " name=\"$txtfile_name\"\n" . 
        "Content-Disposition: attachment;\n" . " filename=\"$txtfile_name\"\n" . 
        "Content-Transfer-Encoding: base64\n\n" . $data . "\n\n";
    $message .= "--{$mime_boundary}\n";

    return $message;
}

function sendEmail($email_content)
{

    // send
    $to = "files@dabam.dx.am";
    $from = "webpage@dabam.dx.am"; 
    $subject ="New work for you"; 
    $headers = "From: $from";

    // headers for attachment 
    $headers .= "\nMIME-Version: 1.0\n" . "Content-Type: multipart/mixed;\n" . " boundary=\"{$mime_boundary}\""; 

    $ok = mail($to, $subject, $email_content, $headers); 
    return $ok;
}


if ($_FILES["dabamfile"]["error"] > 0) 
  {
  echo "Error: " . $_FILES["file"]["error"] . "<br>";
  echo "You have to attach a file.<br>";
  }
else
  {
  //echo "Upload: " . $_FILES["dabamfile"]["name"] . "<br>";
  //echo "Type: " . $_FILES["dabamfile"]["type"] . "<br>";
  //echo "Size: " . ($_FILES["dabamfile"]["size"] / 1024) . " kB<br>";
  //echo "Stored in: " . $_FILES["dabamfile"]["tmp_name"];
  
  $datafile_name = $_FILES["dabamfile"]["tmp_name"];
  $txt_content = createTxtFilecontent($_POST);
  $text_body = "Hello!\n\nHere some new work for you.\nI am very happy to do this to you.\n\nBest regards,\nDABAM";
  $email_content = createEmailContent($text_body,$datafile_name, $txt_content); 
 
  $ok = sendEmail($email_content);
  
  if ($ok) { 
        echo "<p>mail sent!</p>"; 
    } else { 
        echo "<p>mail could not be sent!</p>"; 
    }  
  
  }

  
//foreach($txt_content as $txt_line){
//    print_r($txt_line. "<br />\n");
//}

?>
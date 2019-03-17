<?php

include 'globals.php';
$api_url = $argv[3];
$username = $argv[4];
$password = $argv[5];
$page_name = str_replace("__47__","/",$argv[1]);
echo $argv;
return false;

try {
	$wiki = new Wikimate($api_url);
	if ($wiki -> login($username, $password))
		;
	else {
		$error = $wiki -> getError();
		echo "Wikimate error on login : " . $error['login'] . "\n";
		echo "Have to wait " . $error['wait'] . " seconds\n";
	}
} catch ( Exception $e ) {
	echo "<b>Wikimate error2</b>: " . $e -> getMessage() . "\n";
}
$textpage = file_get_contents($argv[2]);
// create a new page object
$page = $wiki->getPage($argv[1]);
// check if the page exists or not and output the page
$page->setText($textpage);
?>


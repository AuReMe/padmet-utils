<?php

include 'globals.php';
include 'configuration.php';

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
//get the page we want
$page = $wiki->getPage($argv[1]);
$wikiCode= $page->getText();
//Save the page in a file
$fileWikiPage= fopen($argv[1].".txt",'a');
fputs($fileWikiPage,$wikiCode);
fclose($fileWikiPage);
?>


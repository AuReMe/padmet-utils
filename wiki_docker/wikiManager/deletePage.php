<?php

include 'globals.php';
include 'configuration.php';

echo "Begin\n";


try {
	echo "create wiki connection\n";
	$wiki = new Wikimate($api_url);
	echo "try to connect\n";
	if ($wiki -> login($username, $password))
		echo "login ok\n";
	else {
		echo "error\n";
		$error = $wiki -> getError();
		echo "Wikimate error on login : " . $error['login'] . "\n";
		echo "Have to wait " . $error['wait'] . " seconds\n";
	}
} catch ( Exception $e ) {
	echo "<b>Wikimate error2</b>: " . $e -> getMessage() . "\n";
}

// create a new page object
echo "Delete page\n";
$page = $wiki->getPage($argv[1]);
	if ($page->delete("The page was deleted after updating the wiki.")) {
		echo "page $argv[1] deleted";
	} else {
		echo "Cannot delete page $argv[1]";
	}
?>


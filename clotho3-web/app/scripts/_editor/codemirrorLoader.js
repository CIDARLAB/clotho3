
angular.module('clotho.editor').service('codemirrorLoader', function(Clotho, $q) {

	var rootDir = 'bower_components/codemirror/mode/';
	var langMap = {
		javascript : rootDir + 'javascript/javascript.js',
		css : rootDir + 'css/css.js',
		html : rootDir + 'htmlmixed/htmlmixed.js',
		python : rootDir + 'python/python.js',
		groovy : rootDir + 'groovy/groovy.js',
		java : rootDir + 'clike/clike.js',
		markdown : rootDir + 'markdown/markdown.js'
	};

	function loadLanguage (lang) {
		var url = langMap[angular.lowercase(lang)];
		if (url) {
			return $clotho.extensions.mixin(url);
		} else {
			console.log('language not supported in codemirror: ' + lang);
			return $q.when();
		}
	}

	function loadMain () {
		return $q.all([
			$clotho.extensions.mixin('bower_components/codemirror/lib/codemirror.js'),
			$clotho.extensions.css('bower_components/codemirror/lib/codemirror.css')
		]);
	}

	return {
		loadMain : loadMain,
		loadLanguage : loadLanguage
	}

});
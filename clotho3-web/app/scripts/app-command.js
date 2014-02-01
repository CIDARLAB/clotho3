'use strict';

//Command-bar-only stack
angular.module('clothoRoot', ['clotho.core', 'clotho.commandbar'])
.run(function(Clotho) {
	//sort of init() function with server (assumes WebSocket is set up)
	Clotho.submit("clotho.run('clientSetup', [])");
});

angular.element(document).ready(function() {
	angular.bootstrap(document, ['clothoRoot']);
});
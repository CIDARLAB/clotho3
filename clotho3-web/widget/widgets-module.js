'use strict';

Application.Widgets.config(['$controllerProvider', '$compileProvider', '$filterProvider', '$provide', function($controllerProvider, $compileProvider, $filterProvider, $provide) {
    Application.Widgets.providers = {
        $controllerProvider: $controllerProvider,
        $compileProvider: $compileProvider,
        $filterProvider: $filterProvider,
        $provide: $provide
    };

    Application.Widgets.getQueue = function() {
        return angular.module('clotho.widgets')._invokeQueue;
    };

    Application.Widgets.registeredQueue = Application.Widgets.getQueue().length;

}])
    .run(['$rootScope', function($rootScope) {

        //need to call this before compiling new element
        Application.Widgets.processQueue = function() {
            var queue = Application.Widgets.getQueue();
            console.log("processing");

            for(var i=Application.Widgets.registeredQueue;i<queue.length;i++) {
                var call = queue[i];
                // call is in the form [providerName, providerFunc, providerArguments]
                var provider = Application.Widgets.providers[call[0]];
                if(provider) {
                    // e.g. $controllerProvider.register("Ctrl", function() { ... })
                    console.log("provider exists");
                    provider[call[1]].apply(provider, call[2]);
                }
            }
            Application.Widgets.registeredQueue = Application.Widgets.getQueue().length;
        };

        Application.Widgets.recompile = function(element) {
            if (typeof element == 'undefined') {return;}
            $($clotho.appRoot).injector().invoke(function($compile, $rootScope) {
                $compile($(element))($rootScope);
                $rootScope.$apply();
            });
        };

        Application.Widgets.mixin = function(url, element) {
            $script(url, function() {
                Application.Widgets.processQueue();
                Application.Widgets.recompile(element);
            });
        };
}]);

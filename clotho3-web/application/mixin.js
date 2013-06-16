'use strict';

Application.Extensions.config(['$controllerProvider', '$compileProvider', '$filterProvider', '$provide', function($controllerProvider, $compileProvider, $filterProvider, $provide) {

    Application.Extensions.providers = {
        $controllerProvider: $controllerProvider,
        $compileProvider: $compileProvider,
        $filterProvider: $filterProvider,
        $provide: $provide
    };

    Application.Extensions.getQueue = function() {
        return angular.module('clotho.extensions')._invokeQueue;
    };

    Application.Extensions.registeredQueue = Application.Extensions.getQueue().length;

}])
    .run(['$rootScope', '$q', function($rootScope, $q) {

        //need to call this before compiling new element
        Application.Extensions.processQueue = function() {
            var queue = Application.Extensions.getQueue();
            //console.log("processing");

            for(var i=Application.Extensions.registeredQueue;i<queue.length;i++) {
                var call = queue[i];
                // call is in the form [providerName, providerFunc, providerArguments]
                var provider = Application.Extensions.providers[call[0]];
                if(provider) {
                    // e.g. $controllerProvider.register("Ctrl", function() { ... })
                    //console.log("provider exists");
                    provider[call[1]].apply(provider, call[2]);
                }
            }

            //verify make sure this is passed by reference i.e. not only changed locally
            Application.Extensions.registeredQueue = Application.Extensions.getQueue().length;
        };

        Application.Extensions.recompile = function(element) {
            if (typeof element == 'undefined') {return;}
            $($clotho.appRoot).injector().invoke(function($compile, $rootScope) {
                $compile($(element))($rootScope);
                $rootScope.$apply();
            });
        };

        Application.mixin = function(urls, element) {

            var deferred = $q.defer();

            $script(urls, function() {
                Application.Extensions.processQueue();
                Application.Extensions.recompile(element);
                $rootScope.$safeApply(deferred.resolve(element));
            });

            return deferred.promise;
        };

}]);

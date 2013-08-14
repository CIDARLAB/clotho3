'use strict';

//note - ngModel doesn't work for directives with isolate scope..
//  - https://github.com/angular/angular.js/issues/1924
//todo - broaden for other function use
//note - this implementation is meant for directives that already handle ngModel.$render (e.g. input)
Application.Foundation.directive('clothoRun', ['Clotho', '$timeout', function(Clotho, $timeout) {
    return {
        restrict : 'A',
        require : 'ngModel',
        scope: true,
        link: function (scope, element, attrs, ngModel) {

            console.log(scope);
            console.log(ngModel);

            //command, args
            scope.$watch(function() {
                return attrs.clothoRun
            }, function(newval, oldval) {
                console.log(newval);


            });

            //model listeners

            scope.$watch(
            //'model'
            function() {return ngModel.$modelValue}
            , function(newval, oldval) {

                console.log(newval);
                if (!newval) return;

                Clotho.run('blah', newval).then(function(data) {
                    console.log(data);
                    ngModel.$setViewValue(data);
                });
            });


            /*$timeout(function() {
                ngModel.$setViewValue('sup')
            }, 2000);*/


            /*
            var fn = function(input) {

                var promise = Clotho.run('blah', input).then(function (data) { return data });
                console.log(promise);
                return promise;
            };

            ngModel.$formatters.push(fn)
            */

        }
    }
}]);
Application.Trails.directive('trailPage', ['$timeout', '$controller', function($timeout, $controller) {

    // fields that are handled:
    // backend: CSS (url), mixin (array|url), script (array|url), onload (array|url), controller (name, must be mixed in)
    // content: text (html), video (object), template (url), quiz (object), markdown (text)

    return {
        restrict: 'A',
        template: '<div ng-repeat="comp in pageComponents">' +
            '<div trail-page-component="comp"></div>' +
            '</div>',
        scope: {
            page: '=trailPage'
        },
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs) {

                    //future in ng-1.2.x, use notify callbacks for updates
                    //future - this provides the foundation for Clotho.view() -- move it there
                    //todo - error callbacks
                    //todo - flag that component has been mixed in

                    if (!!scope.page.dictionary) {
                        angular.extend(scope, scope.page.dictionary);
                    }

                    return Application.css(scope.page.css)
                        .then(function() {
                            return Application.mixin(scope.page.mixin)
                        })
                        .then(function() {
                            return Application.script(scope.page.script)
                        })
                        .then(function (){

                            // verify this is the best way to do this
                            // todo - create empty as default and use ng-controller in template?
                            //check for controller, must be already included (e.g. by mixin)
                            if (scope.page.controller) {
                                var locals = {};
                                locals.$scope = scope;
                                var ctrl = $controller(scope.page.controller, locals);
                                element.data('$ngControllerController', ctrl);
                            }
                        })
                        .then(function() {

                            /*
                             var promises = [];
                             angular.forEach(scope.page.contents, function (component) {
                                promises.push(loadPageComponent(component.type, component.params))
                             });
                             return promises;
                             */

                            scope.pageComponents = scope.page.contents;
                        }, function(error) {
                            //todo - handle correctly
                            console.log(error);
                        });

                },
                post: function postLink(scope, element, attrs) {
                    if (!!scope.page.onload) {
                        console.log('loading page onload script');
                        $timeout(function() {
                            Application.script(scope.page.onload)
                        });
                    }
                }
            }
        }
    }
}]);
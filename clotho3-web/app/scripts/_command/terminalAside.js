/*
This just serves as a fairly static aside, for use as the terminal

can be activated on Clotho.PubSub channel 'toggleTerminalActive'

a much more extensible version is here:
https://github.com/mgcrea/angular-strap/blob/master/src/aside/aside.js
 */

angular.module('clotho.commandbar')
	.value('terminalAsideOptions', {
		visible : false
	})
	.directive('terminalAside', function($http, $q, $templateCache, $window, $animate, $compile, terminalAsideOptions, hotkeys, Clotho) {

	return {
		restrict: 'EA',
		templateUrl : '../../views/_interface/terminalAside.html',
		replace: true,
		scope: {
			title: '=?asideTitle',
			contentUrl: '=?asideContentUrl'
		},
		controller : function ($scope, $element, $attrs) {
			$scope.submit = function () {
				if ($scope.terminalQuery) {
					Clotho.submit($scope.terminalQuery);
				}
			}
		},
		link: function postLink(scope, element, attrs) {

			function fetchTemplate (url) {
				return $q.when( $templateCache.get(url) || $http.get(url) )
				.then(function (res) {
					if(angular.isObject(res)) {
						$templateCache.put(url, res.data);
						return res.data;
					}
					return res;
				});
			}

			scope.$hide = function () {
				terminalAsideOptions.visible = false;
			};

			scope.$show = function () {
				terminalAsideOptions.visible = true;
			};

			scope.$toggle = function () {
				terminalAsideOptions.visible = !terminalAsideOptions.visible;
				angular.element('body').attr('aside-status', !!terminalAsideOptions.visible ? 'active' : '');
			};

			scope.$watch('contentUrl', function (newval, oldval) {
				if (!!newval) {
					fetchTemplate(newval).then(function (template) {
						console.info('ASIDE TEMPLATE' + template);
						scope.content = $compile(template)(scope);
					});
				}
			});

			Clotho.listen('toggleTerminalActive', scope.$toggle, scope);

			scope.$watch(function () {
				return terminalAsideOptions.visible
			}, function (newval) {
				!!newval ? element.addClass('active') : element.removeClass('active');
			});

			hotkeys.add('t', "Toggle Clotho Terminal", scope.$toggle)

		}
	};
})
.directive('terminalAsideTrigger', function (terminalAsideOptions) {
	return {
		restrict: 'A',
		template: '<div id="terminalAsideTrigger" ng-click="toggle()" ng-attr-status="{{activeClass ? \'active\' : \'\'}}" ng-class="{active : activeClass}"></div>',
		replace: true,
		scope: true,
		link: function postLink(scope, element, attrs, transclusion) {
			scope.toggle = function() {
				terminalAsideOptions.visible = !terminalAsideOptions.visible;
				scope.activeClass = terminalAsideOptions.visible;
				angular.element('body').attr('aside-status', !!terminalAsideOptions.visible ? 'active' : '');
			};

		}
	}
});
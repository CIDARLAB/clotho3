/*
This just serves as a fairly static aside, for use with the terminal
a much more extensible version is here:
https://github.com/mgcrea/angular-strap/blob/master/src/aside/aside.js
 */

angular.module('clotho.interface')
	.value('terminalAsideOptions', {
		visible : false
	})
	.directive('terminalAside', function($http, $q, $templateCache, $window, $animate, $compile, terminalAsideOptions) {

	return {
		restrict: 'EA',
		templateUrl : 'views/_interface/aside.html',
		replace: true,
		scope: {
			title: '=asideTitle',
			contentUrl: '=asideContentUrl'
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

			scope.$watch('contentUrl', function (newval, oldval) {
				if (!!newval) {
					fetchTemplate(newval).then(function (template) {
						console.error('TEMPLATE' + template);
						scope.content = $compile(template)(scope);
					});
				}
			});

			scope.$watch(function () {
				return terminalAsideOptions.visible
			}, function (newval) {
				!!newval ? element.addClass('active') : element.removeClass('active');
			});

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
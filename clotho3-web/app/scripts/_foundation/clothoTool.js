angular.module('clotho.foundation')
.directive('clothoTool', function ($http, $compile, Debug) {

	var Debugger = new Debug('clothoTool', '#cc9999');

	var toolRoot = 'partials/tools/';

	//fixme -- this is so ghetto
	var tools = {
		"revcomp" : {
			partial : "revcomp.html"
		},
		"digestCuts" : {
			partial : "digestCuts.html",
			dependencies : {
				mixin : toolRoot + "digestCuts.js"
			}
		}
	};

	return {
		restrict : "A",
		replace : true,
		link : function (scope, element, attrs) {
			if (angular.isUndefined(attrs.clothoTool)) {
				Debugger.warn('clothoTool attr is not defined');
				return;
			}

			element.addClass('clotho-tool loading');

			scope.$watch(function () {
				return attrs.clothoTool;
			}, function (newval) {
				if (tools[newval]) {
					constructTool(newval);
				} else {
					Debugger.warn('tool name not provided');
				}
			});


			function constructTool (name) {

				$clotho.extensions.downloadDependencies(tools[name].dependencies)
				.then(function (onloadFn) {
					$http.get(toolRoot + tools[name].partial)
					.success(function (data) {
						element.removeClass('loading');
						element.html($compile(data)(scope));
						onloadFn();
					})
					.error(function (data) {
						Debugger.error('Could not find tool' + attrs.clothoTool);
						element.remove();
					});
				});
			}

		}
	}
});
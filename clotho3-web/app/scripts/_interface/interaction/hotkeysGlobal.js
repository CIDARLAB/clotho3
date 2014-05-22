angular.module('clotho.interface')
	.config(function (hotkeysProvider) {
		hotkeysProvider.template = '' +
			'<div class="modal fade hotkeyModal" ng-class="{in: helpVisible}">' +
				'<div class="modal-dialog">' +
					'<div class="modal-content">' +
						'<div class="modal-header">' +
							'<button type="button" class="close" ng-click="closeCheatSheet()">&times;</button>' +
							'<h4 class="modal-title">{{ title }}</h4>' +
						'</div>' +
						'<div class="modal-body">' +
							'<table><tbody>' +
								'<tr ng-repeat="hotkey in hotkeys | filter:{ description: \'!$$undefined$$\' }">' +
									'<td class="hotkey-keys">' +
									'<span ng-repeat="key in hotkey.format() track by $index" class="hotkey-key">{{ key }}</span>' +
									'</td>' +
									'<td class="hotkey-text">{{ hotkey.description }}</td>' +
								'</tr>' +
							'</tbody></table>' +
						'</div>' +
					'</div><!-- /.modal-content -->' +
				'</div><!-- /.modal-dialog -->' +
			'</div><!-- /.modal -->';
	})
.run(function (hotkeys, $location, CommandBar) {
	hotkeys.add('f', 'Focus Command Bar', function (event) {
		event.preventDefault();
		CommandBar.focusInput();
	});
	hotkeys.add('a', 'Show Activity Log', function (event) {
		event.preventDefault();
		CommandBar.display.toggleActivityLog();
	});
	hotkeys.add('g h', 'Go to Homepage', function () {
		$location.path('/')
	});
	hotkeys.add('g b', 'Go to Browser', function () {
		$location.path('/browser')
	});
	hotkeys.add('g e', 'Go to Editor', function () {
		$location.path('/editor')
	});
	hotkeys.add('g t', 'Go to Trails', function () {
		$location.path('/trails')
	});
});
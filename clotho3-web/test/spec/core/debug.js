describe('clotho.core Debug', function() {

	var Debug,
		$log,
		Debugger;

	beforeEach(function() {
		module('clotho.core');

		inject(function (_Debug_, _$log_) {
			Debug = _Debug_;
			$log = _$log_;
		});

		Debugger = new Debug('Module Name', '#ff0000');
	});

	it('should return a factory', function () {
		expect(angular.isObject(Debugger)).toBe(true);
	});

	it('should expose angular $log service', function () {
		expect(Debugger.$log).toEqual($log);
	});

	it('should call $log for some functions - log, warn, error, debug, info', function () {
		spyOn($log, 'log').and.callThrough();
		Debugger.log('something');
		expect($log.log).toHaveBeenCalled();
	});
});